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

#include <algorithm>
#include <cstring>
#include <ctime>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

namespace {

std::string quadruples_record_name(const std::string& label, int ijkl, int component = -1) {
    std::string name = label + " " + std::to_string(ijkl);
    if (component >= 0) name += " " + std::to_string(component);
    return name;
}

int quadruplet_key(int i, int j, int k, int l, int nocc) {
    return ((i * nocc + j) * nocc + k) * nocc + l;
}

void save_quads_record(const std::shared_ptr<PSIO>& psio, const std::string& name, size_t rows, size_t cols,
                       const double* data) {
    const size_t size = rows * cols * sizeof(double);
    // PSIO's legacy write interface takes a mutable buffer but does not modify it.
#pragma omp critical(dlpno_quads_psio)
    psio->write_entry(PSIF_DLPNO_QUADS, name, reinterpret_cast<char*>(const_cast<double*>(data)), size);
}

void load_quads_record(const std::shared_ptr<PSIO>& psio, const std::string& name, size_t rows, size_t cols,
                       double* data) {
    const size_t size = rows * cols * sizeof(double);
#pragma omp critical(dlpno_quads_psio)
    psio->read_entry(PSIF_DLPNO_QUADS, name, reinterpret_cast<char*>(data), size);
}

}  // namespace

DLPNOCCSDT_Q::DLPNOCCSDT_Q(SharedWavefunction ref_wfn, Options &options)
    : DLPNOCCSDT(ref_wfn, options),
      write_quadruples_intermediates_(options.get_bool("WRITE_QUADRUPLES_INTERMEDIATES")) {}
DLPNOCCSDT_Q::~DLPNOCCSDT_Q() {}

void DLPNOCCSDT_Q::print_header() {
    const double t_cut_qno = options_.get_double("T_CUT_QNO");

    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                  DLPNO-CCSDT(Q)              \n");
    outfile->Printf("                   by Andy Jiang               \n");
    outfile->Printf("              DOI: 10.1063/5.0257672          \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n",
                    options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Inherited CCSD and CCSDT parameters are reported in the preceding headers.\n");
    outfile->Printf("  Quadruples parameters:\n");
    outfile->Printf("    T_CUT_QNO                       = %6.3e \n", t_cut_qno);
    outfile->Printf("    T_CUT_QNO_DIAG_SCALE            = %6.3e \n",
                    options_.get_double("T_CUT_QNO_DIAG_SCALE"));
    outfile->Printf("    T_CUT_QNO_STRONG                = %6.3e \n",
                    t_cut_qno * options_.get_double("T_CUT_QNO_STRONG_SCALE"));
    outfile->Printf("    T_CUT_QNO_WEAK                  = %6.3e \n",
                    t_cut_qno * options_.get_double("T_CUT_QNO_WEAK_SCALE"));
    outfile->Printf("    T_CUT_QNO_PRE                   = %6.3e \n",
                    options_.get_double("T_CUT_QNO_PRE"));
    outfile->Printf("    T_CUT_QUADS_WEAK                = %6.3e \n",
                    options_.get_double("T_CUT_QUADS_WEAK"));
    outfile->Printf("    T_CUT_DO_QUADS                  = %6.3e \n",
                    options_.get_double("T_CUT_DO_QUADS"));
    outfile->Printf("    T_CUT_DO_QUADS_PRE              = %6.3e \n",
                    options_.get_double("T_CUT_DO_QUADS_PRE"));
    outfile->Printf("    T_CUT_MKN_QUADS                 = %6.3e \n",
                    options_.get_double("T_CUT_MKN_QUADS"));
    outfile->Printf("    T_CUT_MKN_QUADS_PRE             = %6.3e \n",
                    options_.get_double("T_CUT_MKN_QUADS_PRE"));
    outfile->Printf("    F_CUT_Q                         = %6.3e \n",
                    options_.get_double("F_CUT_Q"));
    outfile->Printf("    T_CUT_ITER_Q                    = %6.3e \n",
                    options_.get_double("T_CUT_ITER_Q"));
    outfile->Printf("    MIN_QNOS                        = %6d   \n",
                    options_.get_int("MIN_QNOS"));
    outfile->Printf("    QUADS_MAX_WEAK_PAIRS            = %6d   \n",
                    options_.get_int("QUADS_MAX_WEAK_PAIRS"));
    outfile->Printf("    E_CONVERGENCE                   = %6.3e \n",
                    options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE                   = %6.3e \n",
                    options_.get_double("R_CONVERGENCE"));
    outfile->Printf("    DLPNO_MAXITER                   = %6d   \n",
                    options_.get_int("DLPNO_MAXITER"));
    outfile->Printf("    Q0_APPROXIMATION                = %6s   \n",
                    options_.get_bool("Q0_APPROXIMATION") ? "TRUE" : "FALSE");
    outfile->Printf("    DLPNO_TOGGLE_MEMORY             = %6s   \n",
                    toggle_memory_ ? "TRUE" : "FALSE");
    outfile->Printf("    WRITE_QUADRUPLES_INTERMEDIATES  = %6s   \n\n",
                    write_quadruples_intermediates_ ? "TRUE" : "FALSE");
}

void DLPNOCCSDT_Q::save_quadruples_tensor(const Tensor<double, 4>& tensor, const std::string& label,
                                          int ijkl) {
    const size_t nqno = static_cast<size_t>(n_qno_[ijkl]);
    const size_t nqno2 = nqno * nqno;
    save_quads_record(psio_, quadruples_record_name(label, ijkl), nqno2, nqno2, tensor.data());
}

Tensor<double, 4> DLPNOCCSDT_Q::load_quadruples_tensor(const std::string& label, int ijkl) {
    const size_t nqno = static_cast<size_t>(n_qno_[ijkl]);
    const size_t nqno2 = nqno * nqno;
    Tensor<double, 4> tensor(label, nqno, nqno, nqno, nqno);
    load_quads_record(psio_, quadruples_record_name(label, ijkl), nqno2, nqno2, tensor.data());
    return tensor;
}

void DLPNOCCSDT_Q::save_quadruplet_energy_intermediates(
    const QuadrupletEnergyIntermediates& intermediates, int ijkl) {
    const size_t nqno = static_cast<size_t>(n_qno_[ijkl]);
    const size_t nlmo = lmoquadruplet_to_lmos_[ijkl].size();
    for (size_t component = 0; component < intermediates.K_iabe.size(); ++component) {
        save_quads_record(psio_, quadruples_record_name("K_iabe", ijkl, component), nqno, nqno * nqno,
                          intermediates.K_iabe[component].data());
    }
    for (size_t component = 0; component < intermediates.K_iajm.size(); ++component) {
        save_quads_record(psio_, quadruples_record_name("K_iajm", ijkl, component), nqno, nlmo,
                          intermediates.K_iajm[component].data());
        save_quads_record(psio_, quadruples_record_name("K_iajb", ijkl, component), nqno, nqno,
                          intermediates.K_iajb[component].data());
        save_quads_record(psio_, quadruples_record_name("U_iajb", ijkl, component), nqno, nqno,
                          intermediates.U_iajb[component].data());
    }
}

DLPNOCCSDT_Q::QuadrupletEnergyIntermediates DLPNOCCSDT_Q::load_quadruplet_energy_intermediates(
    int ijkl) {
    const size_t nqno = static_cast<size_t>(n_qno_[ijkl]);
    const size_t nlmo = lmoquadruplet_to_lmos_[ijkl].size();
    QuadrupletEnergyIntermediates intermediates;
    for (size_t component = 0; component < intermediates.K_iabe.size(); ++component) {
        intermediates.K_iabe[component] = Tensor<double, 3>("K_iabe", nqno, nqno, nqno);
        load_quads_record(psio_, quadruples_record_name("K_iabe", ijkl, component), nqno, nqno * nqno,
                          intermediates.K_iabe[component].data());
    }
    for (size_t component = 0; component < intermediates.K_iajm.size(); ++component) {
        intermediates.K_iajm[component] = Tensor<double, 2>("K_iajm", nqno, nlmo);
        intermediates.K_iajb[component] = Tensor<double, 2>("K_iajb", nqno, nqno);
        intermediates.U_iajb[component] = Tensor<double, 2>("U_iajb", nqno, nqno);
        load_quads_record(psio_, quadruples_record_name("K_iajm", ijkl, component), nqno, nlmo,
                          intermediates.K_iajm[component].data());
        load_quads_record(psio_, quadruples_record_name("K_iajb", ijkl, component), nqno, nqno,
                          intermediates.K_iajb[component].data());
        load_quads_record(psio_, quadruples_record_name("U_iajb", ijkl, component), nqno, nqno,
                          intermediates.U_iajb[component].data());
    }
    return intermediates;
}

Tensor<double, 4> DLPNOCCSDT_Q::matmul_4d(const Tensor<double, 4> &A, const SharedMatrix &X, 
                    int dim_old, int dim_new, bool contract_first) {
    /* Transform all four indices of a rank-4 tensor: A'[i,j,k,l] = A[I,J,K,L] X[i,I] X[j,J] X[k,K] X[l,L]. */

    // TODO: Change this into a TensorView
    Tensor<double, 2> Xview("Xview", dim_new, dim_old);
    ::memcpy(Xview.data(), X->get_pointer(), dim_new * dim_old * sizeof(double));

    Tensor<double, 4> A_new2("A_new2", dim_new, dim_old, dim_old, dim_new);

    if (contract_first) {
        Tensor<double, 4> A_new1("A_new1", dim_old, dim_old, dim_old, dim_new);
        einsum(0.0, Indices{index::I, index::J, index::K, index::l}, &A_new1, 1.0, Indices{index::I, index::J, index::K, index::L}, A, Indices{index::l, index::L}, Xview);
        einsum(0.0, Indices{index::i, index::J, index::K, index::l}, &A_new2, 1.0, Indices{index::I, index::J, index::K, index::l}, A_new1, Indices{index::i, index::I}, Xview);
    } else {
        einsum(0.0, Indices{index::i, index::J, index::K, index::l}, &A_new2, 1.0, Indices{index::i, index::J, index::K, index::L}, A, Indices{index::l, index::L}, Xview);
    }

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
    // Direct occupied/QNO permutation map P_(i,j,k,l), manuscript Eqs. (51)-(55).
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

    int max_weak_pairs = options_.get_int("QUADS_MAX_WEAK_PAIRS");

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

                if (weak_pair_count > max_weak_pairs) continue;

                for (int l : lmopair_to_lmos_[ij]) {
                    int kl = i_j_to_ij_[k][l];
                    if (kl == -1) continue;
                    if (i > l || j > l || k > l) continue;
                    int il_weak = i_j_to_ij_weak_[i][l], jl_weak = i_j_to_ij_weak_[j][l], kl_weak = i_j_to_ij_weak_[k][l];
                    if (i == j && j == l || j == k && k == l || i == k && k == l) continue;

                    if (il_weak != -1) weak_pair_count += 1;
                    if (jl_weak != -1) weak_pair_count += 1;
                    if (kl_weak != -1) weak_pair_count += 1;

                    if (weak_pair_count > max_weak_pairs) continue;

                    ijkl_to_i_j_k_l_.push_back(std::make_tuple(i, j, k, l));

                    std::vector<int> ijkl_list = {i, j, k, l};
                    for (const auto &perm : quadruple_permutations_) {
                        auto &[i_idx, j_idx, k_idx, l_idx] = perm;
                        int ip = ijkl_list[i_idx], jp = ijkl_list[j_idx], kp = ijkl_list[k_idx], lp = ijkl_list[l_idx];
                        i_j_k_l_to_ijkl_[quadruplet_key(ip, jp, kp, lp, naocc)] = ijkl;
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
                for (const auto &perm : quadruple_permutations_) {
                    auto &[i_idx, j_idx, k_idx, l_idx] = perm;
                    int ip = ijkl_list[i_idx], jp = ijkl_list[j_idx], kp = ijkl_list[k_idx], lp = ijkl_list[l_idx];
                    i_j_k_l_to_ijkl_new[quadruplet_key(ip, jp, kp, lp, naocc)] = ijkl_new;
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
        // lmoquadruplet_to_paos_ext_.clear();
    }

    lmoquadruplet_to_ribfs_.resize(n_lmo_quadruplets);
    lmoquadruplet_to_lmos_.resize(n_lmo_quadruplets);
    lmoquadruplet_to_paos_.resize(n_lmo_quadruplets);
    // lmoquadruplet_to_paos_ext_.resize(n_lmo_quadruplets);

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

        // => Make extended domain <= //

        /*
        lmoquadruplet_to_paos_ext_[ijkl] = lmoquadruplet_to_paos_[ijkl];

        for (int m_ijkl = 0; m_ijkl < lmoquadruplet_to_lmos_[ijkl].size(); ++m_ijkl) {
            int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];

            for (int n_ijkl = 0; n_ijkl < lmoquadruplet_to_lmos_[ijkl].size(); ++n_ijkl) {
                int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];

                lmoquadruplet_to_paos_ext_[ijkl] = merge_lists(lmoquadruplet_to_paos_ext_[ijkl], merge_lists(lmo_to_paos[m], lmo_to_paos[n]));

            } // end n_ijkl
        } // end m_ijkl
        */
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
                int ijkl_dense = quadruplet_key(i, j, k, l, naocc);
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

    double t_cut_qno_diag_scale = options_.get_double("T_CUT_QNO_DIAG_SCALE");

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

        // Construct the six-pair quadruplet density of manuscript Eq. (41).
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
        if (i == j || i == k || i == l || j == k || j == l || k == l) qno_scale *= t_cut_qno_diag_scale;

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

    // Allowed unique occupied-index patterns are (1+1+1+1), (2+2), and
    // (2+1+1): C(N,4) + C(N,2) + 3 C(N,3). Triply repeated indices are excluded.
    size_t n_total_possible = (naocc) * (naocc - 1) * (naocc - 2) * (naocc - 3) / 24 + (naocc) * (naocc - 1) / 2 
        + 3 * (naocc) * (naocc - 1) * (naocc - 2) / 6;

    outfile->Printf("  \n");
    outfile->Printf("    Number of (Unique) Local MO quadruplets: %d\n", n_lmo_quadruplets);
    outfile->Printf("    Max Number of Possible (Unique) LMO Quadruplets: %zu (Ratio: %.4f)\n", n_total_possible,
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
    const size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    const size_t n_lmo_pairs = ij_to_i_j_.size();
    const size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    const size_t nthreads = static_cast<size_t>(nthread_);

    size_t qno_basis_memory = 0;
    size_t energy_intermediate_memory = 0;
    size_t quadruples_tensor_memory = 0;
    size_t gamma_workspace_per_thread = 0;
    size_t iteration_workspace_per_thread = 0;

    size_t max_nqno = 0;
    for (size_t ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        max_nqno = std::max(max_nqno, static_cast<size_t>(n_qno_[ijkl]));
    }
    const size_t max_nqno2 = max_nqno * max_nqno;
    const size_t max_nqno4 = max_nqno2 * max_nqno2;

    for (size_t ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        const size_t naux = lmoquadruplet_to_ribfs_[ijkl].size();
        const size_t nlmo = lmoquadruplet_to_lmos_[ijkl].size();
        const size_t npao = lmoquadruplet_to_paos_[ijkl].size();
        const size_t nqno = n_qno_[ijkl];
        const size_t nqno2 = nqno * nqno;
        const size_t nqno3 = nqno2 * nqno;
        const size_t nqno4 = nqno2 * nqno2;
        const size_t nlmo2 = nlmo * nlmo;

        // PAO-to-QNO transformations and QNO orbital energies, Eq. (41).
        qno_basis_memory += npao * nqno + nqno;

        // Algorithm 1 quantities reused by Algorithms 3 and 4 in every energy evaluation.
        const size_t energy_intermediates =
            4 * nqno3 + 16 * nqno * nlmo + 32 * nqno2;
        energy_intermediate_memory += energy_intermediates;
        quadruples_tensor_memory += 2 * nqno4;  // Gamma and T4

        // Conservative peak for one Algorithm 2 task. The estimate includes the sliced
        // DF tensors of Eqs. (47)-(50), all six canonical-Eq. (19) contraction families,
        // and the Algorithm 3/4 energy workspace that is live before the task returns.
        const size_t df_workspace =
            8 * naux * nlmo + 4 * naux * (npao + nqno) +
            2 * naux * nlmo * nqno + 2 * naux * nqno2 + 2 * naux * naux +
            npao * npao;
        const size_t gamma_contraction_workspace =
            16 * nlmo * nqno3 + 16 * nlmo2 + 12 * nlmo * nqno2 +
            32 * nqno2 + 16 * naux * nqno2 + 5 * nqno4;
        const size_t energy_contraction_workspace =
            16 * nlmo * nqno3 + 7 * nqno4 + 8 * nlmo * nqno + 8 * nqno3;
        gamma_workspace_per_thread =
            std::max(gamma_workspace_per_thread,
                     df_workspace + energy_intermediates + gamma_contraction_workspace +
                         energy_contraction_workspace);

        // Algorithm 5 holds the current residual, one or more transformed neighboring T4
        // tensors, and (during the energy pass) the Algorithm 3/4 work arrays. If the
        // persistent bundle is disk-backed, one quadruplet's bundle is loaded per thread.
        const size_t residual_workspace = 8 * std::max(nqno4, max_nqno4);
        const size_t disk_load_workspace =
            write_quadruples_intermediates_ ? energy_intermediates + 2 * nqno4 : 0;
        iteration_workspace_per_thread =
            std::max(iteration_workspace_per_thread,
                     std::max(residual_workspace,
                              energy_contraction_workspace + disk_load_workspace));
    }

    // Standalone CCSDT(Q) releases the CCSD/CCSDT residual intermediates before this stage.
    // Count only the pair/triplet amplitudes and orbital transforms still needed to build
    // QNOs, Gamma, and the energy. CCSDTQ, in contrast, must retain the complete lower-rank
    // state for its subsequent coupled T1/T2/T3/T4 iterations.
    size_t retained_lower_rank_memory = 0;
    if (algorithm_ == DLPNOMethod::CCSDTQ) {
        retained_lower_rank_memory = ccsdt_resident_memory_doubles_;
    } else {
        retained_lower_rank_memory = qij_memory_ + qia_memory_ + qab_memory_;

        for (size_t ij = 0; ij < n_lmo_pairs; ++ij) {
            const size_t npao = lmopair_to_paos_[ij].size();
            const size_t npno = n_pno_[ij];
            retained_lower_rank_memory += npao * npno + npno + 2 * npno * npno;
        }
        for (size_t ijk = 0; ijk < n_lmo_triplets; ++ijk) {
            const size_t npao = lmotriplet_to_paos_[ijk].size();
            const size_t ntno = n_tno_[ijk];
            retained_lower_rank_memory += npao * ntno + ntno + ntno * ntno * ntno;
        }
    }

    auto memory_peaks = [&]() {
        const size_t disk_eligible_resident =
            write_quadruples_intermediates_ ? 0 : energy_intermediate_memory + quadruples_tensor_memory;
        const size_t resident_memory =
            retained_lower_rank_memory + qno_basis_memory + disk_eligible_resident;
        const size_t build_peak = resident_memory + nthreads * gamma_workspace_per_thread;
        const size_t iteration_peak = resident_memory + nthreads * iteration_workspace_per_thread;
        return std::make_tuple(resident_memory, build_peak, iteration_peak);
    };

    const double DOUBLES_TO_GB = 1.0e-9 * sizeof(double);
    const double BYTES_TO_GB = 1.0e-9;

    auto print_estimate = [&](const char* title) {
        const auto [resident_memory, build_peak, iteration_peak] = memory_peaks();
        auto print_memory_line = [&](const std::string& label, size_t words) {
            outfile->Printf("    %-48s : %8.3f [GB]\n", label.c_str(), words * DOUBLES_TO_GB);
        };

        outfile->Printf("\n  ==> %s <==\n\n", title);
        print_memory_line(
            algorithm_ == DLPNOMethod::CCSDTQ
                ? "Retained DLPNO-CCSDT state (for CCSDTQ)"
                : "Retained lower-rank state (standalone CCSDT(Q))",
            retained_lower_rank_memory);
        print_memory_line("QNO transforms and orbital energies", qno_basis_memory);
        print_memory_line("Reusable Algorithm 1 energy intermediates",
                          write_quadruples_intermediates_ ? 0 : energy_intermediate_memory);
        print_memory_line("In-core Gamma and T4 tensors",
                          write_quadruples_intermediates_ ? 0 : quadruples_tensor_memory);
        print_memory_line("Algorithm 2 workspace per thread (" + std::to_string(nthreads) + ")",
                          gamma_workspace_per_thread);
        print_memory_line("Algorithm 5 workspace per thread (" + std::to_string(nthreads) + ")",
                          iteration_workspace_per_thread);
        print_memory_line("Estimated resident memory", resident_memory);
        print_memory_line("Estimated Gamma/(Q0) build peak", build_peak);
        print_memory_line("Estimated iterative-(Q) peak", iteration_peak);
        outfile->Printf("    %-48s : %8.3f [GB]\n\n", "Total memory given", memory_ * BYTES_TO_GB);
    };

    print_estimate("DLPNO-CCSDT(Q) Memory Requirements");

    auto [resident_memory, build_peak, iteration_peak] = memory_peaks();
    size_t required_memory = std::max(build_peak, iteration_peak);
    if (toggle_memory_ && !write_quadruples_intermediates_ &&
        required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total required memory is more than 90%% of available memory.\n");
        outfile->Printf(
            "    Switching Gamma, T4, and reusable iterative-(Q) energy intermediates to disk...\n");
        write_quadruples_intermediates_ = true;

        // The disk-backed iteration workspace includes one loaded bundle per active thread.
        iteration_workspace_per_thread = 0;
        for (size_t ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
            const size_t nlmo = lmoquadruplet_to_lmos_[ijkl].size();
            const size_t nqno = n_qno_[ijkl];
            const size_t nqno2 = nqno * nqno;
            const size_t nqno3 = nqno2 * nqno;
            const size_t nqno4 = nqno2 * nqno2;
            const size_t energy_intermediates =
                4 * nqno3 + 16 * nqno * nlmo + 32 * nqno2;
            const size_t energy_contraction_workspace =
                16 * nlmo * nqno3 + 7 * nqno4 + 8 * nlmo * nqno + 8 * nqno3;
            const size_t residual_workspace = 8 * std::max(nqno4, max_nqno4);
            iteration_workspace_per_thread =
                std::max(iteration_workspace_per_thread,
                         std::max(residual_workspace + 2 * nqno4,
                                  energy_contraction_workspace + energy_intermediates + nqno4));
        }

        std::tie(resident_memory, build_peak, iteration_peak) = memory_peaks();
        required_memory = std::max(build_peak, iteration_peak);
        print_estimate("Updated DLPNO-CCSDT(Q) Memory Requirements");
    }

    if (toggle_memory_ && required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf(
            "  Total required memory remains more than 90%% of available memory after all safe toggles.\n");
        throw PSIEXCEPTION("Too little memory given for the DLPNO-CCSDT(Q) algorithm.");
    }

    if (write_quadruples_intermediates_) {
        outfile->Printf(
            "    Writing Gamma, T4, and reusable iterative-(Q) energy intermediates to disk.\n");
        if (algorithm_ == DLPNOMethod::CCSDTQ) {
            outfile->Printf(
                "    T4 amplitudes will be restored to memory before the full CCSDTQ stage.\n");
        }
        outfile->Printf("\n");
    } else {
        outfile->Printf("    Keeping all persistent iterative-(Q) quantities in core.\n\n");
    }
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

    // The Algorithm 1 energy tensors are owned by one bundle per quadruplet.
    // For (Q0), each bundle dies with its task. Iterative (Q) retains it either
    // in core or in PSIF_DLPNO_QUADS according to estimate_memory().
    quadruplet_energy_intermediates_.clear();
    e_ijkl_.assign(n_lmo_quadruplets, 0.0);

    gamma_ijkl_.clear();
    T_iajbkcld_.clear();
    if (store_amplitudes && !write_quadruples_intermediates_) {
        quadruplet_energy_intermediates_.resize(n_lmo_quadruplets);
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

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Algorithm 1: form the local DF quantities of manuscript Eqs. (47)-(50).

        // number of LMOs in the quadruplet domain
        const int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        // number of PAOs in the quadruplet domain (before removing linear dependencies)
        const int npao_ijkl = lmoquadruplet_to_paos_[ijkl].size();
        // number of auxiliary functions in the quadruplet domain
        const int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        // number of quadruplet natural orbitals in quadruplet domain
        const int nqno_ijkl = n_qno_[ijkl];

        /// => Necessary integrals <= ///

        // (Q_ijkl | i m_ijkl), Eq. (47)
        auto q_io = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);
        auto q_jo = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);
        auto q_ko = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);
        auto q_lo = std::make_shared<Matrix>(naux_ijkl, nlmo_ijkl);

        // (Q_ijkl | i a_ijkl), Eq. (48)
        auto q_iv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);
        auto q_jv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);
        auto q_kv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);
        auto q_lv = std::make_shared<Matrix>(naux_ijkl, npao_ijkl);

        // (Q_ijkl | m_ijkl e_ijkl) and (Q_ijkl | a_ijkl b_ijkl), Eqs. (49)-(50)
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

        constexpr int n_occupied_positions = 4;
        QuadrupletEnergyIntermediates energy_intermediates;

        // Algorithm 1 intermediates for canonical Eq. (19), term 1
        for (int idx = 0; idx < n_occupied_positions; ++idx) {
            energy_intermediates.K_iabe[idx] =
                Tensor<double, 3>("(ia|be)", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::a, index::b, index::e},
                   &energy_intermediates.K_iabe[idx], 1.0,
                   Indices{index::Q, index::a}, q_iv_list[idx],
                   Indices{index::Q, index::b, index::e}, q_vv_ein);
        }

        // Algorithm 1 intermediates for canonical Eq. (19), term 2
        std::array<Tensor<double, 4>, 16> T_mkl_list;
        for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
            const int k = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
                const int l = ijkl_list[j_idx];
                
                // Form K_iajm intermediate
                const int ij_position = i_idx * n_occupied_positions + j_idx;
                energy_intermediates.K_iajm[ij_position] =
                    Tensor<double, 2>("(ia|jm)", nqno_ijkl, nlmo_ijkl);
                einsum(0.0, Indices{index::a, index::m},
                       &energy_intermediates.K_iajm[ij_position], 1.0,
                       Indices{index::Q, index::a}, q_iv_list[i_idx],
                       Indices{index::Q, index::m}, q_io_list[j_idx]);

                // Form T_mkl intermediate
                T_mkl_list[ij_position] =
                    Tensor<double, 4>("T_mkl", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                T_mkl_list[ij_position].zero();

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    int mkl_dense = m * naocc * naocc + k * naocc + l;
                    if (i_j_k_to_ijk_.count(mkl_dense)) {
                        int mkl = i_j_k_to_ijk_[mkl_dense];
                        auto S_ijkl_mkl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[mkl]);
                        S_ijkl_mkl = linalg::triplet(X_qno_[ijkl], S_ijkl_mkl, X_tno_[mkl], true, false, false);
                        auto T_mkl = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[mkl], m, k, l), 
                                                        S_ijkl_mkl, n_tno_[mkl], n_qno_[ijkl]);

                        ::memcpy(&T_mkl_list[ij_position](m_ijkl, 0, 0, 0), T_mkl.data(),
                                 nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                    } // end if
                } // end m_ijkl
            } // end j_idx
        } // end i_idx

        // Algorithm 1 intermediates for canonical Eq. (19), term 3
        std::array<Tensor<double, 2>, 16> K_minj_list;
        for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
            for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
                const int ij_position = i_idx * n_occupied_positions + j_idx;
                K_minj_list[ij_position] =
                    Tensor<double, 2>("(mi|nj)", nlmo_ijkl, nlmo_ijkl);
                einsum(0.0, Indices{index::m, index::n}, &K_minj_list[ij_position], 1.0,
                       Indices{index::Q, index::m}, q_io_list[i_idx],
                       Indices{index::Q, index::n}, q_io_list[j_idx]);
            }
        }

        std::array<Tensor<double, 3>, 4> T_mkac_list;
        for (int idx = 0; idx < n_occupied_positions; ++idx) {
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

        // Algorithm 1 intermediates for canonical Eq. (19), term 4
        std::array<Tensor<double, 3>, 4> K_iame_list;
        for (int idx = 0; idx < n_occupied_positions; ++idx) {
            K_iame_list[idx] = Tensor<double, 3>("(ia|me)", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::a, index::m, index::e}, &K_iame_list[idx], 1.0, Indices{index::Q, index::a}, q_iv_list[idx], Indices{index::Q, index::m, index::e}, q_ov_ein);
        }

        std::array<Tensor<double, 2>, 16> T_ijab_list;
        for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j];
                // if (i_idx > j_idx) continue;

                T_ijab_list[i_idx * n_occupied_positions + j_idx] = Tensor<double, 2>("T_ijab", nqno_ijkl, nqno_ijkl);

                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);
                auto T_ij = linalg::triplet(S_ijkl_ij, T_iajb_[ij], S_ijkl_ij, false, false, true);
                ::memcpy(T_ijab_list[i_idx * n_occupied_positions + j_idx].data(),
                         T_ij->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            } // end j_idx
        } // end i_idx

        // Algorithm 1 intermediates for canonical Eq. (19), term 5
        std::array<Tensor<double, 3>, 4> K_mibe_list;
        for (int idx = 0; idx < n_occupied_positions; ++idx) {
            K_mibe_list[idx] = Tensor<double, 3>("(mi|be)", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::m, index::b, index::e}, &K_mibe_list[idx], 1.0, Indices{index::Q, index::m}, q_io_list[idx], Indices{index::Q, index::b, index::e}, q_vv_ein);
        }

        // Algorithm 1 theta intermediate for canonical Eq. (19), term 6
        std::array<Tensor<double, 3>, 16> theta_Qab;
        for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
            for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
                const int ij_idx = i_idx * n_occupied_positions + j_idx;

                theta_Qab[ij_idx] = Tensor<double, 3>("theta_Qab", naux_ijkl, nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::Q, index::a, index::b}, &theta_Qab[ij_idx], 1.0, Indices{index::Q, index::a, index::e}, q_vv_ein, Indices{index::e, index::b}, T_ijab_list[ij_idx]);
            }
        }

        // Algorithm 1 intermediates reused by the [Q] and (Q) energies, Eqs. (25)-(26)
        for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j], ij_idx = i_idx * n_occupied_positions + j_idx;

                // K_iajb
                energy_intermediates.K_iajb[ij_idx] = Tensor<double, 2>("K_iajb", nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::a, index::b},
                       &energy_intermediates.K_iajb[ij_idx], 1.0,
                       Indices{index::Q, index::a}, q_iv_list[i_idx],
                       Indices{index::Q, index::b}, q_iv_list[j_idx]);

                // U_iajb
                energy_intermediates.U_iajb[ij_idx] = Tensor<double, 2>("U_iajb", nqno_ijkl, nqno_ijkl);
                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);
                auto U_ij_psi = linalg::triplet(S_ijkl_ij, Tt_iajb_[ij], S_ijkl_ij, false, false, true);
                ::memcpy(energy_intermediates.U_iajb[ij_idx].data(),
                         U_ij_psi->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            }
        }

        // => Form gamma_ijkl for all unique (i, j, k, l) permutations in quadruplet ijkl
        Tensor<double, 4> gamma_ijkl("gamma_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
        gamma_ijkl.zero();

        // Algorithm 2 reuses the same contraction when repeated occupied indices make
        // several positional permutations equivalent. Apply all corresponding adjoint
        // permutations immediately so only one rank-four source tensor is resident.
        std::unordered_set<int> computed_permutations;
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            const int occupied_permutation = quadruplet_key(i, j, k, l, naocc);
            if (!computed_permutations.insert(occupied_permutation).second) return;

            const int ij_idx = i_idx * n_occupied_positions + j_idx;
            const int kl_idx = k_idx * n_occupied_positions + l_idx;
            const int kj_idx = k_idx * n_occupied_positions + j_idx;

            Tensor<double, 4> gamma_ijkl_perm("gamma_ijkl_perm", nqno_ijkl, nqno_ijkl,
                                               nqno_ijkl, nqno_ijkl);
            gamma_ijkl_perm.zero();
            Tensor<double, 4> gamma_ijkl_buffer_a("gamma_ijkl_buffer_a", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 4> gamma_ijkl_buffer_b("gamma_ijkl_buffer_b", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);

            // Canonical Eq. (19), term 1; DLPNO Algorithm 2:
            // (i'a|be) t_{j'k'l'}^{ecd}.
            const int jkl_dense = j * naocc * naocc + k * naocc + l;
            if (i_j_k_to_ijk_.count(jkl_dense)) {
                const int jkl = i_j_k_to_ijk_[jkl_dense];
                auto S_ijkl_jkl =
                    submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl],
                                            lmotriplet_to_paos_[jkl]);
                S_ijkl_jkl =
                    linalg::triplet(X_qno_[ijkl], S_ijkl_jkl, X_tno_[jkl], true, false, false);
                auto T_jkl =
                    matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[jkl], j, k, l),
                                      S_ijkl_jkl, n_tno_[jkl], n_qno_[ijkl]);
                einsum(0.0, Indices{index::a, index::b, index::c, index::d},
                       &gamma_ijkl_perm, 1.0, Indices{index::a, index::b, index::e},
                       energy_intermediates.K_iabe[i_idx],
                       Indices{index::e, index::c, index::d}, T_jkl);
            }

            // Canonical Eq. (19), term 2; DLPNO Algorithm 2:
            // -(i'a|j'm) t_{mk'l'}^{bcd}.
            einsum(1.0, Indices{index::a, index::b, index::c, index::d},
                   &gamma_ijkl_perm, -1.0, Indices{index::a, index::m},
                   energy_intermediates.K_iajm[ij_idx],
                   Indices{index::m, index::b, index::c, index::d}, T_mkl_list[kl_idx]);

            // Canonical Eq. (19), term 3; DLPNO Algorithm 2:
            // +(mi'|nj') t_{mk'}^{ac} t_{nl'}^{bd}.
            Tensor<double, 3> gamma_term3("gamma_term3", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::n, index::a, index::c}, &gamma_term3, 1.0,
                   Indices{index::m, index::n}, K_minj_list[ij_idx],
                   Indices{index::m, index::a, index::c}, T_mkac_list[k_idx]);
            einsum(0.0, Indices{index::a, index::c, index::b, index::d},
                   &gamma_ijkl_buffer_a, 1.0, Indices{index::n, index::a, index::c},
                   gamma_term3, Indices{index::n, index::b, index::d}, T_mkac_list[l_idx]);
            permute(Indices{index::a, index::b, index::c, index::d},
                    &gamma_ijkl_buffer_b, Indices{index::a, index::c, index::b, index::d},
                    gamma_ijkl_buffer_a);
            gamma_ijkl_perm += gamma_ijkl_buffer_b;

            // Canonical Eq. (19), term 4; DLPNO Algorithm 2:
            // -2(i'a|me) t_{k'j'}^{eb} t_{ml'}^{cd}.
            Tensor<double, 3> gamma_term4("gamma_term4", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::a, index::m, index::b}, &gamma_term4, 1.0,
                   Indices{index::a, index::m, index::e}, K_iame_list[i_idx],
                   Indices{index::e, index::b}, T_ijab_list[kj_idx]);
            Tensor<double, 3> gamma_term4_transposed("gamma_term4_transposed", nqno_ijkl,
                                                     nqno_ijkl, nlmo_ijkl);
            permute(Indices{index::a, index::b, index::m}, &gamma_term4_transposed,
                    Indices{index::a, index::m, index::b}, gamma_term4);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d},
                   &gamma_ijkl_perm, -2.0, Indices{index::a, index::b, index::m},
                   gamma_term4_transposed, Indices{index::m, index::c, index::d},
                   T_mkac_list[l_idx]);

            // Canonical Eq. (19), term 5; DLPNO Algorithm 2:
            // -2(be|mi') t_{k'j'}^{ce} t_{ml'}^{ad}.
            Tensor<double, 3> gamma_term5("gamma_term5", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::m, index::b, index::c}, &gamma_term5, 1.0,
                   Indices{index::m, index::b, index::e}, K_mibe_list[i_idx],
                   Indices{index::c, index::e}, T_ijab_list[kj_idx]);
            einsum(0.0, Indices{index::a, index::d, index::b, index::c},
                   &gamma_ijkl_buffer_a, 1.0, Indices{index::m, index::a, index::d},
                   T_mkac_list[l_idx], Indices{index::m, index::b, index::c}, gamma_term5);
            permute(Indices{index::a, index::b, index::c, index::d},
                    &gamma_ijkl_buffer_b, Indices{index::a, index::d, index::b, index::c},
                    gamma_ijkl_buffer_a);
            gamma_ijkl_buffer_b *= 2.0;
            gamma_ijkl_perm -= gamma_ijkl_buffer_b;

            // Canonical Eq. (19), term 6; DLPNO Algorithm 2:
            // +(cf|ae) t_{i'j'}^{eb} t_{k'l'}^{fd}.
            einsum(1.0, Indices{index::a, index::b, index::c, index::d},
                   &gamma_ijkl_perm, 1.0, Indices{index::Q, index::a, index::b},
                   theta_Qab[ij_idx], Indices{index::Q, index::c, index::d},
                   theta_Qab[kl_idx]);

            einsums::for_sequence<24UL>([&](auto target_perm_idx) {
                auto &[target_i_idx, target_j_idx, target_k_idx, target_l_idx] =
                    quadruple_permutations_[target_perm_idx];
                const int target_permutation =
                    quadruplet_key(ijkl_list[target_i_idx], ijkl_list[target_j_idx],
                                    ijkl_list[target_k_idx], ijkl_list[target_l_idx], naocc);
                if (target_permutation != occupied_permutation) return;

                permute(Indices{index::a, index::b, index::c, index::d},
                        &gamma_ijkl_buffer_a, std::get<target_perm_idx>(einsum_indices),
                        gamma_ijkl_perm);
                gamma_ijkl += gamma_ijkl_buffer_a;
            });
        });

        gamma_ijkl *= 0.5;

        // Initial semicanonical T4 amplitudes, canonical Eq. (20) and Algorithm 2.
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
        double e_quad = compute_quadruplet_energy(ijkl, T_ijkl, energy_intermediates);
        e_ijkl_[ijkl] = e_quad;
        E_Q0 += e_quad;

        if (store_amplitudes) {
            if (write_quadruples_intermediates_) {
                save_quadruples_tensor(gamma_ijkl, "Gamma", ijkl);
                save_quadruples_tensor(T_ijkl, "T4", ijkl);
                save_quadruplet_energy_intermediates(energy_intermediates, ijkl);
            } else {
                gamma_ijkl_[ijkl] = std::move(gamma_ijkl);
                T_iajbkcld_[ijkl] = std::move(T_ijkl);
                quadruplet_energy_intermediates_[ijkl] = std::move(energy_intermediates);
            }
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

double DLPNOCCSDT_Q::compute_quadruplet_energy(
    int ijkl, const Tensor<double, 4>& T4,
    const QuadrupletEnergyIntermediates& intermediates) {

    int naocc = i_j_to_ij_.size();

    auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
    std::vector<int> ijkl_list = {i, j, k, l};
    constexpr int n_occupied_positions = 4;

    // number of LMOs in the quadruplet domain
    const int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
    // number of quadruplet natural orbitals in quadruplet domain
    const int nqno_ijkl = n_qno_[ijkl];

    double quadruplet_energy = 0.0;
    std::unordered_map<int, double> e_perm_energy;

    std::array<Tensor<double, 4>, 16> T_ijm_list;
    for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
        int i = ijkl_list[i_idx];
        for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
            int j = ijkl_list[j_idx];

            const int ij_position = i_idx * n_occupied_positions + j_idx;
            T_ijm_list[ij_position] =
                Tensor<double, 4>("T_ijm", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            T_ijm_list[ij_position].zero();

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                int ijm_dense = i * naocc * naocc + j * naocc + m;
                if (i_j_k_to_ijk_.count(ijm_dense)) {
                    int ijm = i_j_k_to_ijk_[ijm_dense];
                    auto S_ijkl_ijm = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[ijm]);
                    S_ijkl_ijm = linalg::triplet(X_qno_[ijkl], S_ijkl_ijm, X_tno_[ijm], true, false, false);
                    auto T_ijm = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[ijm], i, j, m), 
                                                    S_ijkl_ijm, n_tno_[ijm], n_qno_[ijkl]);

                    ::memcpy(&T_ijm_list[ij_position](m_ijkl, 0, 0, 0), T_ijm.data(),
                             nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                } // end if
            } // end m_ijkl
        } // end j_idx
    } // end i_idx
        
    einsums::for_sequence<24UL>([&](auto perm_idx) {
        auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
        int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
        int ijkl_idx = quadruplet_key(i, j, k, l, naocc);

        if (!e_perm_energy.count(ijkl_idx)) {
            // Set up e_perm_energy
            e_perm_energy[ijkl_idx] = 0.0;

            // Get quadruples amplitude
            Tensor<double, 4> T_ijkl = quadruples_permuter(T4, i, j, k, l);

            // Fifth-order [Q] contribution: canonical Eq. (25), DLPNO Algorithm 3.
            // u_{kl}^{ab}K_{ij}^{cd} - 2u_{kl}^{bd}L_{ij}^{ac} + u_{kl}^{cd}L_{ij}^{ab}
            Tensor<double, 2> U_kl = intermediates.U_iajb[n_occupied_positions * k_idx + l_idx];
            Tensor<double, 2> K_ij = intermediates.K_iajb[n_occupied_positions * i_idx + j_idx];
            Tensor<double, 2> L_ij = intermediates.K_iajb[n_occupied_positions * i_idx + j_idx];
            L_ij *= 2.0;
            L_ij -= intermediates.K_iajb[n_occupied_positions * j_idx + i_idx];

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

            // Reusable rank-four permutation workspace for Eqs. (27)-(28).
            Tensor<double, 4> T_buffer("T_buffer", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);

            // Canonical Eq. (27), used by the sixth-order (Q) energy in Eq. (26)
            // and DLPNO Algorithm 4.
            // bar{t}_{ijkl}^{abcd} = -2t_{ijkl}^{abcd} - t_{ijkl}^{cdab} + t_{ijkl}^{bacd}
            Tensor<double, 4> T_bar = T_ijkl;
            T_bar *= -2.0;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::c, index::d, index::a, index::b}, T_ijkl);
            T_bar -= T_buffer;
            permute(Indices{index::a, index::b, index::c, index::d}, &T_buffer, Indices{index::b, index::a, index::c, index::d}, T_ijkl);
            T_bar += T_buffer;

            // Canonical Eq. (28), used by Eq. (26) and DLPNO Algorithm 4.
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

            // Alpha and beta contractions from canonical Eqs. (29)-(30), evaluated
            // in the factorized local form of DLPNO Algorithm 4.

            // => 2 - P_{cd} contributions

            Tensor<double, 4> T_bar_dc("T_bar_dc", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::a, index::b, index::d, index::c}, &T_bar_dc, Indices{index::a, index::b, index::c, index::d}, T_bar);
            Tensor<double, 4> T_tilde_dc("T_tilde_dc", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::a, index::b, index::d, index::c}, &T_tilde_dc, Indices{index::a, index::b, index::c, index::d}, T_tilde);

            Tensor<double, 2> alpha_ijm_buffer("alpha_ijm_buffer", nlmo_ijkl, nqno_ijkl);
            const int ij_position = i_idx * n_occupied_positions + j_idx;
            einsum(0.0, Indices{index::m, index::d}, &alpha_ijm_buffer, 2.0,
                    Indices{index::m, index::a, index::b, index::c}, T_ijm_list[ij_position],
                    Indices{index::a, index::b, index::c, index::d}, T_bar);
            einsum(1.0, Indices{index::m, index::d}, &alpha_ijm_buffer, -1.0,
                    Indices{index::m, index::a, index::b, index::c}, T_ijm_list[ij_position],
                    Indices{index::a, index::b, index::c, index::d}, T_bar_dc);

            Tensor<double, 2> beta_ijm_buffer("beta_ijm_buffer", nlmo_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::m, index::d}, &beta_ijm_buffer, 2.0,
                    Indices{index::m, index::a, index::b, index::c}, T_ijm_list[ij_position],
                    Indices{index::a, index::b, index::c, index::d}, T_tilde);
            einsum(1.0, Indices{index::m, index::d}, &beta_ijm_buffer, -1.0,
                    Indices{index::m, index::a, index::b, index::c}, T_ijm_list[ij_position],
                    Indices{index::a, index::b, index::c, index::d}, T_tilde_dc);

            Tensor<double, 2> K_lk_T("K_lk_T", nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::d}, &K_lk_T, Indices{index::d, index::m},
                    intermediates.K_iajm[l_idx * n_occupied_positions + k_idx]);
            Tensor<double, 2> K_kl_T("K_kl_T", nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::d}, &K_kl_T, Indices{index::d, index::m},
                    intermediates.K_iajm[k_idx * n_occupied_positions + l_idx]);

            e_perm_energy[ijkl_idx] += 2.0 * (linear_algebra::dot(alpha_ijm_buffer, K_lk_T) + linear_algebra::dot(beta_ijm_buffer, K_kl_T));

            // 2 - P_{kl} contributions
            int k_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), k) - lmoquadruplet_to_lmos_[ijkl].begin();
            Tensor<double, 3> T_ijk = T_ijm_list[i_idx * n_occupied_positions + j_idx](k_ijkl, All, All, All);
            int l_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), l) - lmoquadruplet_to_lmos_[ijkl].begin();
            Tensor<double, 3> T_ijl = T_ijm_list[i_idx * n_occupied_positions + j_idx](l_ijkl, All, All, All);

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
            permute(Indices{index::e, index::c, index::d}, &g_kdce_T,
                    Indices{index::d, index::c, index::e}, intermediates.K_iabe[k_idx]);
            Tensor<double, 3> g_ldce_T("g_ldce_T", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::e, index::c, index::d}, &g_ldce_T,
                    Indices{index::d, index::c, index::e}, intermediates.K_iabe[l_idx]);

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
    double e_curr = 0.0, e_prev = 0.0;
    bool e_converged = false, r_converged = false;

    double f_cut = options_.get_double("F_CUT_Q");
    double t_cut_iter = options_.get_double("T_CUT_ITER_Q");

    std::vector<double> e_ijkl_old(n_lmo_quadruplets, 0.0);

    while (!(e_converged && r_converged)) {
        // RMS residual for each LMO quadruplet, used to assess convergence.
        std::vector<double> R_iajbkcld_rms(n_lmo_quadruplets, 0.0);

        std::time_t time_start = std::time(nullptr);

#pragma omp parallel for schedule(dynamic, 1)
        for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
            int ijkl = sorted_quadruplets_[ijkl_sorted];
            auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

            int nqno_ijkl = n_qno_[ijkl];
            if (std::fabs(e_ijkl_[ijkl] - e_ijkl_old[ijkl]) < std::fabs(e_ijkl_old[ijkl] * t_cut_iter)) continue;

            // S integrals
            std::vector<int> quadruplet_ext_domain;
            for (int m = 0; m < naocc; ++m) {
                int ijkm_dense = quadruplet_key(i, j, k, m, naocc);
                int ijml_dense = quadruplet_key(i, j, m, l, naocc);
                int imkl_dense = quadruplet_key(i, m, k, l, naocc);
                int mjkl_dense = quadruplet_key(m, j, k, l, naocc);

                if (l != m && i_j_k_l_to_ijkl_.count(ijkm_dense) && std::fabs((*F_lmo_)(l, m)) >= f_cut) {
                    int ijkm = i_j_k_l_to_ijkl_[ijkm_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkm]);
                }

                if (k != m && i_j_k_l_to_ijkl_.count(ijml_dense) && std::fabs((*F_lmo_)(k, m)) >= f_cut) {
                    int ijml = i_j_k_l_to_ijkl_[ijml_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijml]);
                }

                if (j != m && i_j_k_l_to_ijkl_.count(imkl_dense) && std::fabs((*F_lmo_)(j, m)) >= f_cut) {
                    int imkl = i_j_k_l_to_ijkl_[imkl_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[imkl]);
                }

                if (i != m && i_j_k_l_to_ijkl_.count(mjkl_dense) && std::fabs((*F_lmo_)(i, m)) >= f_cut) {
                    int mjkl = i_j_k_l_to_ijkl_[mjkl_dense];
                    quadruplet_ext_domain = merge_lists(quadruplet_ext_domain, lmoquadruplet_to_paos_[mjkl]);
                }
            }
            auto S_ijkl = submatrix_rows_and_cols(*S_pao_, quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkl]);
            S_ijkl = linalg::doublet(S_ijkl, X_qno_[ijkl], false, false);

            // Algorithm 5 implements canonical Eqs. (23)-(24). Keep only the current
            // source and amplitude in a thread when the persistent tensors are on disk.
            Tensor<double, 4> gamma_ijkl = write_quadruples_intermediates_
                                                   ? load_quadruples_tensor("Gamma", ijkl)
                                                   : gamma_ijkl_[ijkl];
            Tensor<double, 4> T_ijkl = write_quadruples_intermediates_
                                               ? load_quadruples_tensor("T4", ijkl)
                                               : T_iajbkcld_[ijkl];
            auto get_t4 = [&](int coupled_ijkl) {
                return write_quadruples_intermediates_
                           ? load_quadruples_tensor("T4", coupled_ijkl)
                           : T_iajbkcld_[coupled_ijkl];
            };

            Tensor<double, 4> R_ijkl("R_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            R_ijkl = gamma_ijkl;

            for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
                for (int b_ijkl = 0; b_ijkl < nqno_ijkl; ++b_ijkl) {
                    for (int c_ijkl = 0; c_ijkl < nqno_ijkl; ++c_ijkl) {
                        for (int d_ijkl = 0; d_ijkl < nqno_ijkl; ++d_ijkl) {
                            (R_ijkl)(a_ijkl, b_ijkl, c_ijkl, d_ijkl) += T_ijkl(a_ijkl, b_ijkl, c_ijkl, d_ijkl) *
                                ((*e_qno_[ijkl])(a_ijkl) + (*e_qno_[ijkl])(b_ijkl) + (*e_qno_[ijkl])(c_ijkl) + (*e_qno_[ijkl])(d_ijkl) 
                                    - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k) - (*F_lmo_)(l, l));
                        } // end d_ijkl
                    } // end c_ijkl
                } // end b_ijkl
            } // end a_ijkl

            for (int m = 0; m < naocc; ++m) {
                int ijkm_dense = quadruplet_key(i, j, k, m, naocc);
                int ijml_dense = quadruplet_key(i, j, m, l, naocc);
                int imkl_dense = quadruplet_key(i, m, k, l, naocc);
                int mjkl_dense = quadruplet_key(m, j, k, l, naocc);

                if (l != m && i_j_k_l_to_ijkl_.count(ijkm_dense) && std::fabs((*F_lmo_)(l, m)) >= f_cut) {
                    int ijkm = i_j_k_l_to_ijkl_[ijkm_dense];
                    std::vector<int> ijkm_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkm]);
                    auto S_ijkl_ijkm = linalg::doublet(submatrix_rows(*S_ijkl, ijkm_idx_list), X_qno_[ijkm], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(get_t4(ijkm), i, j, k, m),
                                            S_ijkl_ijkm, n_qno_[ijkm], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(l, m);
                    R_ijkl -= T_temp;
                }

                if (k != m && i_j_k_l_to_ijkl_.count(ijml_dense) && std::fabs((*F_lmo_)(k, m)) >= f_cut) {
                    int ijml = i_j_k_l_to_ijkl_[ijml_dense];
                    std::vector<int> ijml_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijml]);
                    auto S_ijkl_ijml = linalg::doublet(submatrix_rows(*S_ijkl, ijml_idx_list), X_qno_[ijml], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(get_t4(ijml), i, j, m, l),
                                            S_ijkl_ijml, n_qno_[ijml], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(k, m);
                    R_ijkl -= T_temp;
                }

                if (j != m && i_j_k_l_to_ijkl_.count(imkl_dense) && std::fabs((*F_lmo_)(j, m)) >= f_cut) {
                    int imkl = i_j_k_l_to_ijkl_[imkl_dense];
                    std::vector<int> imkl_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[imkl]);
                    auto S_ijkl_imkl = linalg::doublet(submatrix_rows(*S_ijkl, imkl_idx_list), X_qno_[imkl], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(get_t4(imkl), i, m, k, l),
                                            S_ijkl_imkl, n_qno_[imkl], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(j, m);
                    R_ijkl -= T_temp;
                }

                if (i != m && i_j_k_l_to_ijkl_.count(mjkl_dense) && std::fabs((*F_lmo_)(i, m)) >= f_cut) {
                    int mjkl = i_j_k_l_to_ijkl_[mjkl_dense];
                    std::vector<int> mjkl_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[mjkl]);
                    auto S_ijkl_mjkl = linalg::doublet(submatrix_rows(*S_ijkl, mjkl_idx_list), X_qno_[mjkl], true, false);
                    auto T_temp = matmul_4d(quadruples_permuter(get_t4(mjkl), m, j, k, l),
                                            S_ijkl_mjkl, n_qno_[mjkl], n_qno_[ijkl]);
                    T_temp *= (*F_lmo_)(i, m);
                    R_ijkl -= T_temp;
                }
            }

            // => Update T4 Amplitudes <= //
            for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
                for (int b_ijkl = 0; b_ijkl < nqno_ijkl; ++b_ijkl) {
                    for (int c_ijkl = 0; c_ijkl < nqno_ijkl; ++c_ijkl) {
                        for (int d_ijkl = 0; d_ijkl < nqno_ijkl; ++d_ijkl) {
                            T_ijkl(a_ijkl, b_ijkl, c_ijkl, d_ijkl) -= R_ijkl(a_ijkl, b_ijkl, c_ijkl, d_ijkl) /
                                ((*e_qno_[ijkl])(a_ijkl) + (*e_qno_[ijkl])(b_ijkl) + (*e_qno_[ijkl])(c_ijkl) + (*e_qno_[ijkl])(d_ijkl) 
                                    - (*F_lmo_)(i, i) - (*F_lmo_)(j, j) - (*F_lmo_)(k, k) - (*F_lmo_)(l, l));
                        } // end d_ijkl
                    } // end c_ijkl
                } // end b_ijkl
            } // end a_ijkl

            if (write_quadruples_intermediates_) {
                save_quadruples_tensor(T_ijkl, "T4", ijkl);
            } else {
                T_iajbkcld_[ijkl] = std::move(T_ijkl);
            }
            
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
            double e_ijkl = 0.0;
            if (write_quadruples_intermediates_) {
                auto T_ijkl = load_quadruples_tensor("T4", ijkl);
                auto energy_intermediates = load_quadruplet_energy_intermediates(ijkl);
                e_ijkl = compute_quadruplet_energy(ijkl, T_ijkl, energy_intermediates);
            } else {
                e_ijkl = compute_quadruplet_energy(
                    ijkl, T_iajbkcld_[ijkl], quadruplet_energy_intermediates_[ijkl]);
            }
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
    DLPNOCCSDT::compute_energy();

    einsums::profile::initialize();
    print_header();

    // A standalone CCSDT(Q) calculation can release lower-rank integral
    // intermediates that a subsequent full-CCSDTQ iteration would still need.
    if (algorithm_ == DLPNOMethod::CCSDT_Q) {
        // Clear CCSD integrals
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

        // Clear CCSDT integrals
        q_io_.clear();
        q_jo_.clear();
        q_ko_.clear();
        q_iv_.clear();
        q_jv_.clear();
        q_kv_.clear();
        q_ov_.clear();
        q_vv_.clear();
    }

    // Re-create T_iajbkc_clone intermediate
    int n_lmo_triplets = ijk_to_i_j_k_.size();
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];

        T_iajbkc_clone_[ijk] = Tensor<double, 3>("T_ijk", n_tno_[ijk], n_tno_[ijk], n_tno_[ijk]);
        ::memcpy(T_iajbkc_clone_[ijk].data(), T_iajbkc_[ijk]->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
    }

    if (algorithm_ == DLPNOMethod::CCSDT_Q) {
        // The standalone method needs only the Einsums T3 clone from this point onward.
        // Full CCSDTQ retains the lower-rank state for its coupled iterations instead.
        W_iajbkc_.clear();
        V_iajbkc_.clear();
        T_iajbkc_.clear();
        T_n_ijk_.clear();
        U_iajbkc_.clear();
    }

    double t_cut_qno_pre = options_.get_double("T_CUT_QNO_PRE");
    double t_cut_qno = options_.get_double("T_CUT_QNO");

    // Step 1: loose-QNO prescreening and screened-quadruplet correction, Eq. (56).
    outfile->Printf("\n   Starting Quadruplet Prescreening...\n");
    outfile->Printf("     T_CUT_QNO set to %6.3e \n", t_cut_qno_pre);
    outfile->Printf("     T_CUT_DO  set to %6.3e \n", options_.get_double("T_CUT_DO_QUADS_PRE"));
    outfile->Printf("     T_CUT_MKN set to %6.3e \n\n", options_.get_double("T_CUT_MKN_QUADS_PRE"));

    quadruples_sparsity(true);
    qno_transform(t_cut_qno_pre);
    double E_Q0_pre = compute_gamma_ijkl(false);
    outfile->Printf("    (Initial) DLPNO-(Q0) Correlation Energy: %16.12f\n\n", E_Q0_pre);

    // Step 2: semicanonical (Q0) energy for the surviving quadruplets, Eq. (57).
    outfile->Printf("\n   Continuing computation with surviving quadruplets...\n");
    outfile->Printf("     Eliminated all quadruplets with energy less than %6.3e Eh... \n\n", options_.get_double("T_CUT_QUADS_WEAK"));
    quadruples_sparsity(false);
    outfile->Printf("    * Energy Contribution From Screened Quadruplets: %.12f \n\n", de_lccsdt_q_screened_);
    
    outfile->Printf("     T_CUT_QNO (re)set to %6.3e \n", options_.get_double("T_CUT_QNO"));
    outfile->Printf("     T_CUT_DO  (re)set to %6.3e \n", options_.get_double("T_CUT_DO_QUADS"));
    outfile->Printf("     T_CUT_MKN (re)set to %6.3e \n\n", options_.get_double("T_CUT_MKN_QUADS"));

    qno_transform(t_cut_qno);
    double E_Q0 = compute_gamma_ijkl(false);
    e_lccsdt_q_ = E_Q0;
    E_Q_ = E_Q0;
    outfile->Printf("    (Total) DLPNO-(Q0) Correlation Energy:      %16.12f\n", E_Q0 + de_lccsdt_q_screened_);
    outfile->Printf("    * Screened Quadruplets Contribution:        %16.12f\n", de_lccsdt_q_screened_);

    double e_scf = variables_["SCF TOTAL ENERGY"];
    double e_ccsdt_q_corr = E_Q0 + de_lccsdt_q_screened_ + e_lccsdt_ + de_lccsd_t_screened_ + de_weak_ +
                            de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
    double e_ccsdt_q_total = e_scf + e_ccsdt_q_corr;

    outfile->Printf("\n\n  @Total DLPNO-CCSDT(Q0) Energy: %16.12f\n", e_ccsdt_q_total);

    double e_total = e_ccsdt_q_total;

    if (!options_.get_bool("Q0_APPROXIMATION")) {
        // Step 3: iterative non-semicanonical (Q), Algorithm 5 and Eq. (62).
        outfile->Printf("\n\n  ==> Computing Full Iterative (Q) <==\n\n");

        double t_cut_qno_strong_scale = options_.get_double("T_CUT_QNO_STRONG_SCALE");
        double t_cut_qno_weak_scale = options_.get_double("T_CUT_QNO_WEAK_SCALE");
        outfile->Printf("     T_CUT_QNO (re)set to %6.3e for strong quadruplets \n", t_cut_qno * t_cut_qno_strong_scale);
        outfile->Printf("     T_CUT_QNO (re)set to %6.3e for weak quadruplets   \n\n", t_cut_qno * t_cut_qno_weak_scale);

        // Sort quadruplets into "strong" and "weak" quadruplets
        sort_quadruplets(E_Q0);
        qno_transform(t_cut_qno);
        estimate_memory();

        if (write_quadruples_intermediates_) {
            psio_->open(PSIF_DLPNO_QUADS, PSIO_OPEN_NEW);
        }

        double E_Q0_iteration_domains = compute_gamma_ijkl(true);
        double E_Q = lccsdt_q_iterations();
        double dE_Q = E_Q - E_Q0_iteration_domains;
        E_Q_ = E_Q;
        e_lccsdt_q_ = E_Q0 + dE_Q;

        if (write_quadruples_intermediates_) {
            if (algorithm_ == DLPNOMethod::CCSDTQ) {
                // The full-quadruples stage directly consumes T_iajbkcld_. Rehydrate only
                // T4; Gamma and the perturbative energy tensors are no longer needed.
                T_iajbkcld_.resize(ijkl_to_i_j_k_l_.size());
#pragma omp parallel for schedule(dynamic, 1)
                for (int ijkl = 0; ijkl < ijkl_to_i_j_k_l_.size(); ++ijkl) {
                    T_iajbkcld_[ijkl] = load_quadruples_tensor("T4", ijkl);
                }
            }
            psio_->close(PSIF_DLPNO_QUADS, 0);
        }

        gamma_ijkl_.clear();
        quadruplet_energy_intermediates_.clear();
        if (algorithm_ == DLPNOMethod::CCSDT_Q) T_iajbkcld_.clear();

        outfile->Printf("\n");
        outfile->Printf("    DLPNO-(Q0) energy at looser tolerance:      %16.12f\n",
                        E_Q0_iteration_domains);
        outfile->Printf("    DLPNO-(Q)  energy at looser tolerance:      %16.12f\n", E_Q);
        outfile->Printf("    * Net Iterative (Q) contribution:           %16.12f\n\n", dE_Q);

        outfile->Printf("    (Total) DLPNO-(Q) Correlation Energy:       %16.12f\n", E_Q0 + dE_Q + de_lccsdt_q_screened_);
        outfile->Printf("    * DLPNO-(Q0) Contribution:                  %16.12f\n", E_Q0);
        outfile->Printf("    * DLPNO-(Q) Contribution:                   %16.12f\n", dE_Q);
        outfile->Printf("    * Screened Quadruplets Contribution:        %16.12f\n", de_lccsdt_q_screened_);

        // Overall DLPNO-CCSDT(Q) assembly, Eq. (63), including the preceding
        // local-pair/triplet truncation corrections.
        e_ccsdt_q_corr = E_Q0 + dE_Q + de_lccsdt_q_screened_ + e_lccsdt_ + de_lccsd_t_screened_ + de_weak_ +
                         de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
        e_ccsdt_q_total = e_scf + e_ccsdt_q_corr;
        e_total = e_ccsdt_q_total;

        outfile->Printf("\n\n  @Total DLPNO-CCSDT(Q) Energy: %16.12f\n", e_ccsdt_q_total);
    }
    outfile->Printf("  DLPNO-CCSDT(Q) computation complete.\n\n");

    set_scalar_variable("CCSDT(Q) CORRELATION ENERGY", e_ccsdt_q_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsdt_q_corr);
    set_scalar_variable("CCSDT(Q) TOTAL ENERGY", e_ccsdt_q_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsdt_q_total);
    set_scalar_variable("(Q) CORRECTION ENERGY",
                        e_ccsdt_q_total - variables_["CCSDT TOTAL ENERGY"]);

    print_results();

    einsums::profile::finalize();

    timer_off("DLPNO-CCSDT(Q)");

    return e_total;
}

void DLPNOCCSDT_Q::print_results() {
    const int naocc = i_j_to_ij_.size();
    double t1_diagnostic = 0.0;
#pragma omp parallel for reduction(+ : t1_diagnostic)
    for (int i = 0; i < naocc; ++i) {
        t1_diagnostic += T_ia_[i]->vector_dot(T_ia_[i]);
    }
    t1_diagnostic = std::sqrt(t1_diagnostic / (2.0 * naocc));
    outfile->Printf("\n  T1 Diagnostic: %8.8f \n", t1_diagnostic);
    // Lee and Taylor's conventional 0.02 closed-shell threshold was defined from CCSD
    // singles amplitudes (Int. J. Quantum Chem. 36, 199, 1989; DOI: 10.1002/qua.560360824).
    // Retain it here as a conservative warning about the single-reference description.
    constexpr double closed_shell_t1_warning = 0.02;
    if (t1_diagnostic > closed_shell_t1_warning) {
        outfile->Printf("    WARNING: T1 Diagnostic is greater than 0.02; single-reference "
                        "coupled-cluster results may be unreliable.\n");
    }
    set_scalar_variable("CC T1 DIAGNOSTIC", t1_diagnostic);

    const bool q0_approximation = options_.get_bool("Q0_APPROXIMATION");
    const char* quadruples_label = q0_approximation ? "DLPNO-(Q0) Contribution" : "DLPNO-(Q) Contribution";
    const double lower_rank_correlation =
        e_lccsdt_ + de_lccsd_t_screened_ + de_weak_ + de_lmp2_eliminated_ +
        de_dipole_ + de_pno_total_;
    const double quadruples_correlation = e_lccsdt_q_ + de_lccsdt_q_screened_;
    const double total_correlation = lower_rank_correlation + quadruples_correlation;
    const double total_energy = variables_["SCF TOTAL ENERGY"] + total_correlation;
    const double ccsdt_q_minus_ccsdt = total_energy - variables_["CCSDT TOTAL ENERGY"];

    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSDT(Q) Correlation Energy: %16.12f \n", total_correlation);
    outfile->Printf("    LCCSDT Correlation Energy:             %16.12f \n", e_lccsdt_);
    outfile->Printf("    %-38s %16.12f \n", quadruples_label, e_lccsdt_q_);
    outfile->Printf("    Screened Quadruplets Contribution:     %16.12f \n", de_lccsdt_q_screened_);
    outfile->Printf("    Screened Triplets Contribution:        %16.12f \n", de_lccsd_t_screened_);
    outfile->Printf("    Weak Pair Contribution:                %16.12f \n", de_weak_);
    outfile->Printf("    Eliminated Pair MP2 Correction:        %16.12f \n", de_lmp2_eliminated_);
    outfile->Printf("    Dipole Pair Correction:                %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:             %16.12f \n", de_pno_total_);
    outfile->Printf("    CCSDT(Q) - CCSDT Energy:               %16.12f \n", ccsdt_q_minus_ccsdt);
    outfile->Printf("\n\n  @Total DLPNO-CCSDT(Q) Energy: %16.12f \n", total_energy);
}

DLPNOCCSDTQ::DLPNOCCSDTQ(SharedWavefunction ref_wfn, Options& options) : DLPNOCCSDT_Q(ref_wfn, options) {}
DLPNOCCSDTQ::~DLPNOCCSDTQ() {}

void DLPNOCCSDTQ::print_header() {
    double t_cut_qno_full = options_.get_double("T_CUT_QNO_FULL");

    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                   DLPNO-CCSDTQ                \n");
    outfile->Printf("                   by Andy Jiang               \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    T_CUT_QNO_FULL             = %6.3e \n", options_.get_double("T_CUT_QNO_FULL"));
}

void DLPNOCCSDTQ::estimate_memory() {

    outfile->Printf("\n ==> DLPNO-CCSDTQ Memory Estimate <== \n\n");

    size_t naocc = i_j_to_ij_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    size_t q_io_memory = 0;
    size_t q_iv_memory = 0;
    size_t q_ov_memory = 0;
    size_t q_vv_memory = 0;

#pragma omp parallel for reduction(+ : q_io_memory, q_iv_memory, q_ov_memory, q_vv_memory)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
        auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

        int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        int npao_ijkl = lmoquadruplet_to_paos_[ijkl].size();
        int nqno_ijkl = n_qno_[ijkl];

        q_io_memory += 4 * naux_ijkl * nlmo_ijkl;
        q_iv_memory += 4 * naux_ijkl * nqno_ijkl;

        if (!disk_ints_quads_) {
            q_ov_memory += naux_ijkl * nlmo_ijkl * nqno_ijkl;
            q_vv_memory += naux_ijkl * nqno_ijkl * nqno_ijkl;
        }
    } // end ijkl_sorted

    if (!disk_ints_quads_) {
        outfile->Printf("    Keeping all LMO/QNO ERIs in core!\n\n");
    } else {
        outfile->Printf("    Writing expensive (Q_{ijkl} | m_{ijkl} a_{ijkl}) and (Q_{ijkl} | a_{ijkl} b_{ijkl}) integrals to disk!\n\n");
    }
}

void DLPNOCCSDTQ::xpno_transform(double xpno_tolerance) {
    timer_on("XPNO transform");

    int naocc = nalpha_ - nfrzc();
    int nbf = basisset_->nbf();
    int n_lmo_pairs = ij_to_i_j_.size();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    int min_pnos = options_.get_int("MIN_PNOS");
    double xpno_diag_scale = options_.get_double("T_CUT_XPNO_DIAG_SCALE");
    double xpno_occ_tolerance = options_.get_double("T_CUT_TRACE_XPNO");

    lmopair_to_paos_ext_.resize(n_lmo_pairs);
    X_pno_ext_.resize(n_lmo_pairs);
    e_pno_ext_.resize(n_lmo_pairs);
    n_pno_ext_.resize(n_lmo_pairs);

    std::vector<double> occ_xpno(n_lmo_pairs, 0.0);
    std::vector<double> trace_xpno(n_lmo_pairs, 0.0);

#pragma omp parallel for
    for (int kl = 0; kl < n_lmo_pairs; ++kl) {
        auto &[k, l] = ij_to_i_j_[kl];
        int lk = ij_to_ji_[kl];

        if (k > l) continue;

        lmopair_to_paos_ext_[kl] = lmopair_to_paos_[kl];

        for (int mn = 0; mn < n_lmo_pairs; ++mn) {
            auto &[m, n] = ij_to_i_j_[mn];

            // int mnkl_idx = m * std::pow(naocc, 3) + n * std::pow(naocc, 2) + k * naocc + l;
            // int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;
            // if (mnkl == -1) continue;

            int mk = i_j_to_ij_[m][k], ml = i_j_to_ij_[m][l], nk = i_j_to_ij_[n][k], nl = i_j_to_ij_[n][l];
            if (mk == -1 || ml == -1 || nk == -1 || nl == -1) continue;

            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[mn]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[mk]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[ml]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[nk]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[nl]);
        }
        lmopair_to_paos_ext_[lk] = lmopair_to_paos_ext_[kl];

        // number of PAOs in the extended domain of kl
        int npao_ext_kl = lmopair_to_paos_ext_[kl].size();

        //                                               //
        // ==> Canonicalize extended PAOs of pair kl <== //
        //                                               //

        auto S_pao_kl_ext = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_ext_[kl]);
        auto F_pao_kl_ext = submatrix_rows_and_cols(*F_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_ext_[kl]);

        SharedMatrix X_pao_kl_ext;
        SharedVector e_pao_kl_ext;
        std::tie(X_pao_kl_ext, e_pao_kl_ext) = orthocanonicalizer(S_pao_kl_ext, F_pao_kl_ext);

        F_pao_kl_ext = linalg::triplet(X_pao_kl_ext, F_pao_kl_ext, X_pao_kl_ext, true, false, false);

        // number of PAOs in the domain after removing linear dependencies
        int npao_can_kl_ext = X_pao_kl_ext->colspi(0);

        //                                            //
        // ==> Canonical PAOs  to Canonical XPNOs <== //
        //                                            //

        size_t nvir_kl_ext = F_pao_kl_ext->rowspi(0);

        // Compute pair density from amplitudes
        SharedMatrix D_kl_sum = std::make_shared<Matrix>("D_kl_sum", nvir_kl_ext, nvir_kl_ext);
        D_kl_sum->zero();

        // Take the sum of the pair density over quadruplets
        int quad_count = 0;
        int nvir_ijkl_avg = 0;
        for (int mn = 0; mn < n_lmo_pairs; ++mn) {
            auto &[m, n] = ij_to_i_j_[mn];

            // int m_kl = lmopair_to_lmos_dense_[kl][m], n_kl = lmopair_to_lmos_dense_[kl][n];
            // if (m_kl == -1 || n_kl == -1) continue;

            int mk = i_j_to_ij_[m][k], ml = i_j_to_ij_[m][l], nk = i_j_to_ij_[n][k], nl = i_j_to_ij_[n][l];
            if (mk == -1 || ml == -1 || nk == -1 || nl == -1) continue;

            // kl
            auto D_kl = linalg::doublet(Tt_iajb_[kl], T_iajb_[kl], false, true);
            D_kl->add(linalg::doublet(Tt_iajb_[kl], T_iajb_[kl], true, false));
            if (k == l) D_kl->scale(0.5);
            auto S_kl = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_[kl]);
            S_kl = linalg::triplet(X_pao_kl_ext, S_kl, X_pno_[kl], true, false, false);
            D_kl_sum->add(linalg::triplet(S_kl, D_kl, S_kl, false, false, true));

            // mn
            auto D_mn = linalg::doublet(Tt_iajb_[mn], T_iajb_[mn], false, true);
            D_mn->add(linalg::doublet(Tt_iajb_[mn], T_iajb_[mn], true, false));
            if (m == n) D_mn->scale(0.5);
            auto S_mn = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_[mn]);
            S_mn = linalg::triplet(X_pao_kl_ext, S_mn, X_pno_[mn], true, false, false);
            D_kl_sum->add(linalg::triplet(S_mn, D_mn, S_mn, false, false, true));

            // mk
            auto D_mk = linalg::doublet(Tt_iajb_[mk], T_iajb_[mk], false, true);
            D_mk->add(linalg::doublet(Tt_iajb_[mk], T_iajb_[mk], true, false));
            if (m == k) D_mk->scale(0.5);
            auto S_mk = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_[mk]);
            S_mk = linalg::triplet(X_pao_kl_ext, S_mk, X_pno_[mk], true, false, false);
            D_kl_sum->add(linalg::triplet(S_mk, D_mk, S_mk, false, false, true));

            // ml
            auto D_ml = linalg::doublet(Tt_iajb_[ml], T_iajb_[ml], false, true);
            D_ml->add(linalg::doublet(Tt_iajb_[ml], T_iajb_[ml], true, false));
            if (m == l) D_ml->scale(0.5);
            auto S_ml = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_[ml]);
            S_ml = linalg::triplet(X_pao_kl_ext, S_ml, X_pno_[ml], true, false, false);
            D_kl_sum->add(linalg::triplet(S_ml, D_ml, S_ml, false, false, true));

            // nk
            auto D_nk = linalg::doublet(Tt_iajb_[nk], T_iajb_[nk], false, true);
            D_nk->add(linalg::doublet(Tt_iajb_[nk], T_iajb_[nk], true, false));
            if (n == k) D_nk->scale(0.5);
            auto S_nk = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_[nk]);
            S_nk = linalg::triplet(X_pao_kl_ext, S_nk, X_pno_[nk], true, false, false);
            D_kl_sum->add(linalg::triplet(S_nk, D_nk, S_nk, false, false, true));

            // nl
            auto D_nl = linalg::doublet(Tt_iajb_[nl], T_iajb_[nl], false, true);
            D_nl->add(linalg::doublet(Tt_iajb_[nl], T_iajb_[nl], true, false));
            if (n == l) D_nl->scale(0.5);
            auto S_nl = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_[nl]);
            S_nl = linalg::triplet(X_pao_kl_ext, S_nl, X_pno_[nl], true, false, false);
            D_kl_sum->add(linalg::triplet(S_nl, D_nl, S_nl, false, false, true));

            int mnkl_idx = m * std::pow(naocc, 3) + n * std::pow(naocc, 2) + k * naocc + l;
            int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;
            if (mnkl != -1) nvir_ijkl_avg += n_qno_[mnkl];
            
            quad_count++;
        }
        D_kl_sum->scale(1.0 / (6.0 * quad_count));

        // Diagonalization of pair density
        auto X_pno_kl = std::make_shared<Matrix>("eigenvectors", nvir_kl_ext, nvir_kl_ext);
        Vector pno_occ("eigenvalues", nvir_kl_ext);
        D_kl_sum->diagonalize(*X_pno_kl, pno_occ, descending);

        // Compute trace sum
        double occ_total = 0.0;
        for (size_t a = 0; a < nvir_kl_ext; ++a) {
            occ_total += pno_occ.get(a);
        }

        int nvir_kl_final = 0;
        double occ_curr = 0.0;

        double xpno_scale = 1.0;
        if (k == l) xpno_scale = xpno_diag_scale;

        for (size_t a = 0; a < nvir_kl_ext; ++a) {
            if (fabs(pno_occ.get(a)) >= xpno_scale * xpno_tolerance || occ_curr / occ_total < xpno_occ_tolerance || a < min_pnos) {
                occ_curr += pno_occ.get(a);
                nvir_kl_final++;
            } // end if
        } // end a

        nvir_kl_final = std::max(min_pnos, nvir_kl_final);

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(nvir_kl_final);

        // This transformation gives orbitals that are orthonormal but not canonical
        X_pno_kl = X_pno_kl->get_block({zero, X_pno_kl->rowspi()}, {zero, dim_final});
        pno_occ = pno_occ.get_block({zero, dim_final});

        SharedMatrix pno_canon;
        SharedVector e_pno_kl;
        std::tie(pno_canon, e_pno_kl) = canonicalizer(X_pno_kl, F_pao_kl_ext);

        // This transformation gives orbitals that are orthonormal and canonical
        X_pno_kl = linalg::doublet(X_pno_kl, pno_canon, false, false);
        X_pno_kl = linalg::doublet(X_pao_kl_ext, X_pno_kl, false, false);

        X_pno_ext_[kl] = X_pno_kl;
        e_pno_ext_[kl] = e_pno_kl;
        n_pno_ext_[kl] = X_pno_kl->colspi(0);
        occ_xpno[kl] = pno_occ.get(n_pno_ext_[kl] - 1);
        trace_xpno[kl] = occ_curr / occ_total;

        // account for symmetry
        if (k < l) {
            X_pno_ext_[lk] = X_pno_kl;
            e_pno_ext_[lk] = e_pno_kl;
            n_pno_ext_[lk] = X_pno_kl->colspi(0);
            occ_xpno[lk] = occ_xpno[kl];
            trace_xpno[lk] = trace_xpno[kl];
        } // end if (k < l)
    }

    // Print out PNO domain information
    int pno_count_total = 0, pno_count_min = nbf, pno_count_max = 0;
    double occ_number_total = 0.0, occ_number_min = 2.0, occ_number_max = 0.0;
    double trace_total = 0.0, trace_min = 1.0, trace_max = 0.0;
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        pno_count_total += n_pno_ext_[ij];
        pno_count_min = std::min(pno_count_min, n_pno_ext_[ij]);
        pno_count_max = std::max(pno_count_max, n_pno_ext_[ij]);
        occ_number_total += occ_xpno[ij];
        occ_number_min = std::min(occ_number_min, occ_xpno[ij]);
        occ_number_max = std::max(occ_number_max, occ_xpno[ij]);
        trace_total += trace_xpno[ij];
        trace_min = std::min(trace_min, trace_xpno[ij]);
        trace_max = std::max(trace_max, trace_xpno[ij]);
    }

    outfile->Printf("  \n");
    outfile->Printf("    (Extended) Natural Orbitals per Local MO pair:\n");
    outfile->Printf("      Avg: %3d NOs \n", pno_count_total / n_lmo_pairs);
    outfile->Printf("      Min: %3d NOs \n", pno_count_min);
    outfile->Printf("      Max: %3d NOs \n", pno_count_max);
    outfile->Printf("      Avg Occ Number Tol: %.3e \n", occ_number_total / n_lmo_pairs);
    outfile->Printf("      Min Occ Number Tol: %.3e \n", occ_number_min);
    outfile->Printf("      Max Occ Number Tol: %.3e \n", occ_number_max);
    outfile->Printf("      Avg Trace Sum: %.6f \n", trace_total / n_lmo_pairs);
    outfile->Printf("      Min Trace Sum: %.6f \n", trace_min);
    outfile->Printf("      Max Trace Sum: %.6f \n", trace_max);
    outfile->Printf("  \n");

    timer_off("XPNO transform");
}

void DLPNOCCSDTQ::compute_integrals() {

    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // Four-center quantities
    q_io_list_.resize(n_lmo_quadruplets);
    q_iv_list_.resize(n_lmo_quadruplets);
    q_ov_ijkl_.resize(n_lmo_quadruplets);
    q_vv_ijkl_.resize(n_lmo_quadruplets);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
        auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

        int nqno_ijkl = n_qno_[ijkl];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        
        // number of auxiliary functions in the quadruplet domain
        const int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        // number of LMOs in the quadruplet domain
        const int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        // number of PAOs in the quadruplet domain
        const int npao_ijkl = lmoquadruplet_to_paos_[ijkl].size();

        std::stringstream q_ov_name;
        q_ov_name << "(Q_ijkl | m a) " << (ijkl);
        std::stringstream q_vv_name;
        q_vv_name << "(Q_ijkl | a b) " << (ijkl);

        // (Q_{ijkl} | [i, j, k, l] m_{ijkl})
        auto q_io = std::make_shared<Matrix>("(Q_ijkl | m i)", naux_ijkl, nlmo_ijkl);
        auto q_jo = std::make_shared<Matrix>("(Q_ijkl | m j)", naux_ijkl, nlmo_ijkl);
        auto q_ko = std::make_shared<Matrix>("(Q_ijkl | m k)", naux_ijkl, nlmo_ijkl);
        auto q_lo = std::make_shared<Matrix>("(Q_ijkl | m l)", naux_ijkl, nlmo_ijkl);

        // (Q_{ijkl} | [i, j, k, l] a_{ijkl})
        auto q_iv = std::make_shared<Matrix>("(Q_ijkl | i a)", naux_ijkl, npao_ijkl);
        auto q_jv = std::make_shared<Matrix>("(Q_ijkl | j b)", naux_ijkl, npao_ijkl);
        auto q_kv = std::make_shared<Matrix>("(Q_ijkl | k c)", naux_ijkl, npao_ijkl);
        auto q_lv = std::make_shared<Matrix>("(Q_ijkl | l d)", naux_ijkl, npao_ijkl);

        auto q_ov = std::make_shared<Matrix>(q_ov_name.str(), naux_ijkl, nlmo_ijkl * nqno_ijkl);
        auto q_vv = std::make_shared<Matrix>(q_vv_name.str(), naux_ijkl, nqno_ijkl * nqno_ijkl);

        for (int q_ijkl = 0; q_ijkl < naux_ijkl; ++q_ijkl) {
            const int q = lmoquadruplet_to_ribfs_[ijkl][q_ijkl];
            const int centerq = ribasis_->function_to_center(q);

            // Cheaper integrals
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
                } // end v_ijkl
            } // end u_ijkl
            q_vv_tmp = linalg::triplet(X_qno_[ijkl], q_vv_tmp, X_qno_[ijkl], true, false, false);
            ::memcpy(&(*q_vv)(q_ijkl, 0), &(*q_vv_tmp)(0, 0), nqno_ijkl * nqno_ijkl * sizeof(double));
        } // end q_ijkl

        // Contract Intermediates
        q_iv = linalg::doublet(q_iv, X_qno_[ijkl]);
        q_jv = linalg::doublet(q_jv, X_qno_[ijkl]);
        q_kv = linalg::doublet(q_kv, X_qno_[ijkl]);
        q_lv = linalg::doublet(q_lv, X_qno_[ijkl]);
        
        // Multiply by (P|Q)^{-1/2}
        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmoquadruplet_to_ribfs_[ijkl], lmoquadruplet_to_ribfs_[ijkl]);
        A_solve->power(0.5, 1.0e-14);

        C_DGESV_wrapper(A_solve->clone(), q_io);
        C_DGESV_wrapper(A_solve->clone(), q_jo);
        C_DGESV_wrapper(A_solve->clone(), q_ko);
        C_DGESV_wrapper(A_solve->clone(), q_lo);
        C_DGESV_wrapper(A_solve->clone(), q_iv);
        C_DGESV_wrapper(A_solve->clone(), q_jv);
        C_DGESV_wrapper(A_solve->clone(), q_kv);
        C_DGESV_wrapper(A_solve->clone(), q_lv);
        C_DGESV_wrapper(A_solve->clone(), q_ov);
        C_DGESV_wrapper(A_solve->clone(), q_vv);

        q_io_list_[ijkl][0] = Tensor<double, 2>("(Q_ijkl | m i)", naux_ijkl, nlmo_ijkl);
        q_io_list_[ijkl][1] = Tensor<double, 2>("(Q_ijkl | m j)", naux_ijkl, nlmo_ijkl);
        q_io_list_[ijkl][2] = Tensor<double, 2>("(Q_ijkl | m k)", naux_ijkl, nlmo_ijkl);
        q_io_list_[ijkl][3] = Tensor<double, 2>("(Q_ijkl | m l)", naux_ijkl, nlmo_ijkl);

        q_iv_list_[ijkl][0] = Tensor<double, 2>("(Q_ijkl | i a)", naux_ijkl, nqno_ijkl);
        q_iv_list_[ijkl][1] = Tensor<double, 2>("(Q_ijkl | j b)", naux_ijkl, nqno_ijkl);
        q_iv_list_[ijkl][2] = Tensor<double, 2>("(Q_ijkl | k c)", naux_ijkl, nqno_ijkl);
        q_iv_list_[ijkl][3] = Tensor<double, 2>("(Q_ijkl | l d)", naux_ijkl, nqno_ijkl);

        q_ov_ijkl_[ijkl] = Tensor<double, 3>(q_ov->name(), naux_ijkl, nlmo_ijkl, nqno_ijkl);
        q_vv_ijkl_[ijkl] = Tensor<double, 3>(q_vv->name(), naux_ijkl, nqno_ijkl, nqno_ijkl);

        ::memcpy(q_io_list_[ijkl][0].data(), q_io->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));
        ::memcpy(q_io_list_[ijkl][1].data(), q_jo->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));
        ::memcpy(q_io_list_[ijkl][2].data(), q_ko->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));
        ::memcpy(q_io_list_[ijkl][3].data(), q_lo->get_pointer(), naux_ijkl * nlmo_ijkl * sizeof(double));

        ::memcpy(q_iv_list_[ijkl][0].data(), q_iv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_iv_list_[ijkl][1].data(), q_jv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_iv_list_[ijkl][2].data(), q_kv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_iv_list_[ijkl][3].data(), q_lv->get_pointer(), naux_ijkl * nqno_ijkl * sizeof(double));

        ::memcpy(q_ov_ijkl_[ijkl].data(), q_ov->get_pointer(), naux_ijkl * nlmo_ijkl * nqno_ijkl * sizeof(double));
        ::memcpy(q_vv_ijkl_[ijkl].data(), q_vv->get_pointer(), naux_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));

        if (disk_ints_quads_) {
#pragma omp critical
            q_ov->save(psio_.get(), PSIF_DLPNO_QIA_QNO, ::psi::Matrix::SubBlocks);

#pragma omp critical
            q_vv->save(psio_.get(), PSIF_DLPNO_QAB_QNO, ::psi::Matrix::ThreeIndexLowerTriangle);

            q_ov_name << "(Q_ijkl | m a) " << (ijkl);
            q_vv_name << "(Q_ijkl | a b) " << (ijkl);

            q_ov_ijkl_[ijkl] = Tensor<double, 3>(q_ov->name(), 0, 0, 0);
            q_vv_ijkl_[ijkl] = Tensor<double, 3>(q_vv->name(), 0, 0, 0);
        }
    }
}

inline Tensor<double, 3> DLPNOCCSDTQ::QIA_QNO(const int ijkl) {
    if (disk_ints_quads_) {
        int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        int nqno_ijkl = n_qno_[ijkl];

        std::stringstream q_ov_name;
        q_ov_name << "(Q_ijkl | m a) " << (ijkl);

        auto q_ov = std::make_shared<Matrix>(q_ov_name.str(), naux_ijkl, nlmo_ijkl * nqno_ijkl);
#pragma omp critical
        q_ov->load(psio_.get(), PSIF_DLPNO_QIA_QNO, ::psi::Matrix::SubBlocks);

        q_ov_ijkl_[ijkl] = Tensor<double, 3>(q_ov->name(), naux_ijkl, nlmo_ijkl, nqno_ijkl);
        ::memcpy(q_ov_ijkl_[ijkl].data(), q_ov->get_pointer(), naux_ijkl * nlmo_ijkl * nqno_ijkl * sizeof(double));
    }

    return q_ov_ijkl_[ijkl];
}

inline Tensor<double, 3> DLPNOCCSDTQ::QAB_QNO(const int ijkl) {
    if (disk_ints_quads_) {
        int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        int nqno_ijkl = n_qno_[ijkl];

        std::stringstream q_vv_name;
        q_vv_name << "(Q_ijkl | a b) " << (ijkl);

        auto q_vv = std::make_shared<Matrix>(q_vv_name.str(), naux_ijkl, nqno_ijkl * nqno_ijkl);
#pragma omp critical
        q_vv->load(psio_.get(), PSIF_DLPNO_QAB_QNO, ::psi::Matrix::ThreeIndexLowerTriangle);

        q_vv_ijkl_[ijkl] = Tensor<double, 3>(q_vv->name(), naux_ijkl, nqno_ijkl, nqno_ijkl);
        ::memcpy(q_vv_ijkl_[ijkl].data(), q_vv->get_pointer(), naux_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
    }

    return q_vv_ijkl_[ijkl];
}

inline Tensor<double, 4> DLPNOCCSDTQ::alpha_ijkl_helper(const Tensor<double, 4>& T_ijkl) {
    int nqno_ijkl = T_ijkl.dim(0);
    Tensor<double, 4> alpha_ijkl = T_ijkl;
    alpha_ijkl *= 2.0;

    for (int a = 0; a < nqno_ijkl; ++a) {
        for (int b = 0; b < nqno_ijkl; ++b) {
            for (int c = 0; c < nqno_ijkl; ++c) {
                for (int d = 0; d < nqno_ijkl; ++d) {
                    alpha_ijkl(a, b, c, d) -= (T_ijkl(b, a, c, d) + T_ijkl(c, b, a, d) + T_ijkl(d, b, c, a));
                } // end d
            } // end c
        } // end b
    } // end a

    return alpha_ijkl;
}

inline Tensor<double, 4> DLPNOCCSDTQ::beta_ijkl_helper(const Tensor<double, 4>& alpha_ijkl) {
    int nqno_ijkl = alpha_ijkl.dim(0);
    Tensor<double, 4> beta_ijkl = alpha_ijkl;
    beta_ijkl *= 2.0;

    for (int a = 0; a < nqno_ijkl; ++a) {
        for (int b = 0; b < nqno_ijkl; ++b) {
            for (int c = 0; c < nqno_ijkl; ++c) {
                for (int d = 0; d < nqno_ijkl; ++d) {
                    beta_ijkl(a, b, c, d) -= (alpha_ijkl(a, c, b, d) + alpha_ijkl(a, d, c, b));
                } // end d
            } // end c
        } // end b
    } // end a

    return beta_ijkl;
}

Tensor<double, 4> DLPNOCCSDTQ::quadruples_spin_summation(const Tensor<double, 4> &X) {

    Tensor<double, 4> alpha = alpha_ijkl_helper(X);
    Tensor<double, 4> beta = beta_ijkl_helper(alpha);
    Tensor<double, 4> gamma = beta;
    gamma *= 2.0;
    gamma -= quadruples_permuter(beta, 0, 1, 3, 2);

    return gamma;
}

Tensor<double, 4> DLPNOCCSDTQ::quadruples_spin_desummation(const Tensor<double, 4> &X) {
    int i = 0, j = 1, k = 2, l = 3;

    // 7/96 term
    Tensor<double, 4> X1 = X;
    X1 *= 35.0;

    // 1/480 terms
    Tensor<double, 4> X2 = quadruples_permuter(X, j, i, k, l);
    X2 += quadruples_permuter(X, k, j, i, l);
    X2 += quadruples_permuter(X, l, j, k, i);
    X2 += quadruples_permuter(X, i, k, j, l);
    X2 += quadruples_permuter(X, i, l, k, j);
    X2 += quadruples_permuter(X, i, j, l, k);
    X2 *= 1.0;

    // 11/480 terms
    Tensor<double, 4> X3 = quadruples_permuter(X, k, i, l, j);
    X3 += quadruples_permuter(X, l, i, j, k);
    X3 += quadruples_permuter(X, j, l, i, k);
    X3 += quadruples_permuter(X, l, k, i, j);
    X3 += quadruples_permuter(X, j, k, l, i);
    X3 += quadruples_permuter(X, k, l, j, i);
    X3 *= 11.0;

    // 1/32 term
    Tensor<double, 4> X4 = quadruples_permuter(X, j, i, l, k);
    X4 += quadruples_permuter(X, k, l, i, j);
    X4 += quadruples_permuter(X, l, k, j, i);
    X4 *= 15.0;

    X1 += X2;
    X1 += X3;
    X1 += X4;
    X1 *= 1.0 / 480.0;

    return X1;
}

void DLPNOCCSDTQ::compute_R_iajb_quads(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb, 
                                        std::vector<std::vector<SharedMatrix>>& R_iajb_buffer) {
    
    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // Compute residual from CCSDT
    DLPNOCCSDT::compute_R_iajb_triples(R_iajb, Rn_iajb, R_iajb_buffer);

    // Clean buffers
#pragma omp parallel for
    for (int thread = 0; thread < nthread_; ++thread) {
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            R_iajb_buffer[thread][ij]->zero();
        } // end ij
    } // end thread

    // Loop over unique quadruplets
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
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

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        std::unordered_set<int> computed_perms;

        // Loop over all possible quadruplet permutations
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            int ijkl_idx = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;
            int kl = i_j_to_ij_[k][l];

            if (!computed_perms.count(ijkl_idx)) {
                // Get quadruples amplitude
                Tensor<double, 4> T_ijkl = quadruples_permuter(T_iajbkcld_[ijkl], i, j, k, l);

                // Form alpha (Jiang and Matthews Equation 3)
                Tensor<double, 4> alpha = alpha_ijkl_helper(T_ijkl);

                // Form beta from alpha (Jiang and Matthews Equation 4)
                Tensor<double, 4> beta = beta_ijkl_helper(alpha);

                // Jiang and Matthews Equation 5
                // R_{kl}^{cd} += 0.25 * P_{kl}^{cd}[(ia|jb) beta_{ijkl}^{abcd}]
                Tensor<double, 2> K_iajb("K_iajb", nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::a, index::b}, &K_iajb, 1.0, Indices{index::Q, index::a}, q_iv_list_[ijkl][i_idx],
                        Indices{index::Q, index::b}, q_iv_list_[ijkl][j_idx]);

                Tensor<double, 2> R_iajb_cont("R_iajb_cont", nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::c, index::d}, &R_iajb_cont, 0.25, Indices{index::a, index::b, index::c, index::d}, beta,
                        Indices{index::a, index::b}, K_iajb);

                // Form QNO overlap integral (S_ijkl_kl)
                auto S_ijkl_kl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[kl]);
                S_ijkl_kl = linalg::triplet(X_qno_[ijkl], S_ijkl_kl, X_pno_[kl], true, false, false);

                // Copy into a Psi4 Matrix for the subsequent Matrix-based basis transformation.
                auto R_iajb_psi = std::make_shared<Matrix>(nqno_ijkl, nqno_ijkl);
                ::memcpy(R_iajb_psi->get_pointer(), R_iajb_cont.data(), nqno_ijkl * nqno_ijkl * sizeof(double));
                R_iajb_buffer[thread][kl]->add(linalg::triplet(S_ijkl_kl, R_iajb_psi, S_ijkl_kl, true, false, false));

                // Add ijkl_idx to computed perms
                computed_perms.emplace(ijkl_idx);
            }
        });
    }

    // Flush buffers
#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int ji = ij_to_ji_[ij];
        for (int thread = 0; thread < nthread_; ++thread) {
            R_iajb[ij]->add(R_iajb_buffer[thread][ij]);
            R_iajb[ij]->add(R_iajb_buffer[thread][ji]->transpose());
        } // end thread
    } // end ij
}

void DLPNOCCSDTQ::compute_R_iajbkc_quads(std::vector<SharedMatrix>& R_iajbkc) {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // Compute residual from CCSDT
    DLPNOCCSDT::compute_R_iajbkc(R_iajbkc);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];

        int ntno_ijk = n_tno_[ijk];
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();

        auto R_ijk = R_iajbkc[ijk];
        auto T_ijk = T_iajbkc_[ijk];

        // Read integrals when disk-backed storage is enabled.
        if (disk_ints_) {
            q_ov_[ijk] = QIA_TNO(ijk);
            q_vv_[ijk] = QAB_TNO(ijk);
        }

        // => STEP 0: T1-dress DF integrals <= //
        auto i_ijk = std::find(lmotriplet_to_lmos_[ijk].begin(), lmotriplet_to_lmos_[ijk].end(), i) - lmotriplet_to_lmos_[ijk].begin();
        auto j_ijk = std::find(lmotriplet_to_lmos_[ijk].begin(), lmotriplet_to_lmos_[ijk].end(), j) - lmotriplet_to_lmos_[ijk].begin();
        auto k_ijk = std::find(lmotriplet_to_lmos_[ijk].begin(), lmotriplet_to_lmos_[ijk].end(), k) - lmotriplet_to_lmos_[ijk].begin();

        Tensor<double, 1> T_i("T_i", ntno_ijk);
        ::memcpy(T_i.data(), &(T_n_ijk_[ijk])(i_ijk, 0), ntno_ijk * sizeof(double));
        Tensor<double, 1> T_j("T_j", ntno_ijk);
        ::memcpy(T_j.data(), &(T_n_ijk_[ijk])(j_ijk, 0), ntno_ijk * sizeof(double));
        Tensor<double, 1> T_k("T_k", ntno_ijk);
        ::memcpy(T_k.data(), &(T_n_ijk_[ijk])(k_ijk, 0), ntno_ijk * sizeof(double));

        Tensor<double, 2> q_iv_t1 = q_iv_[ijk];
        einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1, -1.0, Indices{index::Q, index::l}, q_io_[ijk], Indices{index::l, index::a}, T_n_ijk_[ijk]);
        einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1, 1.0, Indices{index::Q, index::a, index::b}, q_vv_[ijk], Indices{index::b}, T_i);
        Tensor<double, 2> q_iv_t1_temp("q_iv_t1_temp", naux_ijk, nlmo_ijk);
        einsum(0.0, Indices{index::Q, index::l}, &q_iv_t1_temp, 1.0, Indices{index::Q, index::l, index::b}, q_ov_[ijk], Indices{index::b}, T_i);
        einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1, -1.0, Indices{index::Q, index::l}, q_iv_t1_temp, Indices{index::l, index::a}, T_n_ijk_[ijk]);

        Tensor<double, 2> q_jv_t1 = q_jv_[ijk];
        einsum(1.0, Indices{index::Q, index::a}, &q_jv_t1, -1.0, Indices{index::Q, index::l}, q_jo_[ijk], Indices{index::l, index::a}, T_n_ijk_[ijk]);
        einsum(1.0, Indices{index::Q, index::a}, &q_jv_t1, 1.0, Indices{index::Q, index::a, index::b}, q_vv_[ijk], Indices{index::b}, T_j);
        Tensor<double, 2> q_jv_t1_temp("q_jv_t1_temp", naux_ijk, nlmo_ijk);
        einsum(0.0, Indices{index::Q, index::l}, &q_jv_t1_temp, 1.0, Indices{index::Q, index::l, index::b}, q_ov_[ijk], Indices{index::b}, T_j);
        einsum(1.0, Indices{index::Q, index::a}, &q_jv_t1, -1.0, Indices{index::Q, index::l}, q_jv_t1_temp, Indices{index::l, index::a}, T_n_ijk_[ijk]);

        Tensor<double, 2> q_kv_t1 = q_kv_[ijk];
        einsum(1.0, Indices{index::Q, index::a}, &q_kv_t1, -1.0, Indices{index::Q, index::l}, q_ko_[ijk], Indices{index::l, index::a}, T_n_ijk_[ijk]);
        einsum(1.0, Indices{index::Q, index::a}, &q_kv_t1, 1.0, Indices{index::Q, index::a, index::b}, q_vv_[ijk], Indices{index::b}, T_k);
        Tensor<double, 2> q_kv_t1_temp("q_kv_t1_temp", naux_ijk, nlmo_ijk);
        einsum(0.0, Indices{index::Q, index::l}, &q_kv_t1_temp, 1.0, Indices{index::Q, index::l, index::b}, q_ov_[ijk], Indices{index::b}, T_k);
        einsum(1.0, Indices{index::Q, index::a}, &q_kv_t1, -1.0, Indices{index::Q, index::l}, q_kv_t1_temp, Indices{index::l, index::a}, T_n_ijk_[ijk]);

        // This one is special... the second index is dressed instead of the first (per convention), to increase computational efficiency
        Tensor<double, 3> q_vv_t1 = q_vv_[ijk];
        Tensor<double, 3> q_vo("q_vo", naux_ijk, ntno_ijk, nlmo_ijk);
        permute(Indices{index::Q, index::a, index::l}, &q_vo, Indices{index::Q, index::l, index::a}, q_ov_[ijk]);
        einsum(1.0, Indices{index::Q, index::a, index::b}, &q_vv_t1, -1.0, Indices{index::Q, index::a, index::l}, q_vo, Indices{index::l, index::b}, T_n_ijk_[ijk]);

        Tensor<double, 2> q_io_t1 = q_io_[ijk];
        einsum(1.0, Indices{index::Q, index::l}, &q_io_t1, 1.0, Indices{index::Q, index::l, index::a}, q_ov_[ijk], Indices{index::a}, T_i);
        Tensor<double, 2> q_jo_t1 = q_jo_[ijk];
        einsum(1.0, Indices{index::Q, index::l}, &q_jo_t1, 1.0, Indices{index::Q, index::l, index::a}, q_ov_[ijk], Indices{index::a}, T_j);
        Tensor<double, 2> q_ko_t1 = q_ko_[ijk];
        einsum(1.0, Indices{index::Q, index::l}, &q_ko_t1, 1.0, Indices{index::Q, index::l, index::a}, q_ov_[ijk], Indices{index::a}, T_k);

        std::vector<int> ijk_idx = {i, j, k};
        const int THREE = ijk_idx.size();
        std::vector<Tensor<double, 1>> T_i_list = {T_i, T_j, T_k};
        std::vector<Tensor<double, 2>> q_io_list = {q_io_[ijk], q_jo_[ijk], q_ko_[ijk]};
        std::vector<Tensor<double, 2>> q_io_t1_list = {q_io_t1, q_jo_t1, q_ko_t1};
        std::vector<Tensor<double, 2>> q_iv_t1_list = {q_iv_t1, q_jv_t1, q_kv_t1};

        // => Form Fock Matrix Intermediates <= //

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

        // K_menj integrals
        std::array<Tensor<double, 3>, 3> K_menj_list;
        for (int j_idx = 0; j_idx < THREE; ++j_idx) {
            K_menj_list[j_idx] = Tensor<double, 3>("K_menj", ntno_ijk, nlmo_ijk, nlmo_ijk);
            einsum(0.0, Indices{index::e, index::m, index::n}, &K_menj_list[j_idx], 1.0, Indices{index::Q, index::e, index::m}, q_vo,
                    Indices{index::Q, index::n}, q_io_t1_list[j_idx]);
        }

        // K_dble integrals
        Tensor<double, 4> K_ledb("K_ledb", nlmo_ijk, ntno_ijk, ntno_ijk, ntno_ijk);
        einsum(0.0, Indices{index::l, index::e, index::d, index::b}, &K_ledb, 1.0, Indices{index::Q, index::l, index::e}, q_ov_[ijk], Indices{index::Q, index::d, index::b}, q_vv_t1);

        // Permutations
        std::vector<std::tuple<int, int, int>> short_perms_idx = {std::make_tuple(0, 1, 2), std::make_tuple(1, 0, 2), std::make_tuple(2, 1, 0)};

        // Total residual contribution
        Tensor<double, 3> R_ijk_cont("R_ijk_cont", ntno_ijk, ntno_ijk, ntno_ijk);
        
        for (int perm_idx = 0; perm_idx < short_perms_idx.size(); ++perm_idx) {
            auto &[i_idx, j_idx, k_idx] = short_perms_idx[perm_idx];
            int i = ijk_idx[i_idx], j = ijk_idx[j_idx], k = ijk_idx[k_idx];

            Tensor<double, 3> R_ijk_buffer("R_ijk_buffer", ntno_ijk, ntno_ijk, ntno_ijk);
            R_ijk_buffer.zero();

            for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
                int m = lmotriplet_to_lmos_[ijk][m_ijk];
                int mijk_idx = m * std::pow(naocc, 3) + i * std::pow(naocc, 2) + j * naocc + k;
                int mijk = i_j_k_l_to_ijkl_.count(mijk_idx) ? (i_j_k_l_to_ijkl_[mijk_idx]) : -1;

                if (mijk != -1) {
                    auto S_ijk_mijk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], lmoquadruplet_to_paos_[mijk]);
                    S_ijk_mijk = linalg::triplet(X_tno_[ijk], S_ijk_mijk, X_qno_[mijk], true, false, false);
                    Tensor<double, 2> S_ijk_mijk_ein("S_ijk_mijk_ein", n_tno_[ijk], n_qno_[mijk]);
                    ::memcpy(S_ijk_mijk_ein.data(), S_ijk_mijk->get_pointer(), n_tno_[ijk] * n_qno_[mijk] * sizeof(double));

                    Tensor<double, 4> T_mijk = quadruples_permuter(T_iajbkcld_[mijk], m, i, j, k);
                    Tensor<double, 4> alpha_mijk = alpha_ijkl_helper(T_mijk);

                    // Jiang and Matthews Eq. 6a 1/6 F_{me} alpha_{mijk}^{eabc}
                    Tensor<double, 1> F_me_slice = F_ld(m_ijk, All);
                    Tensor<double, 1> F_me_qno("F_me_qno", n_qno_[mijk]);
                    einsum(0.0, Indices{index::E}, &F_me_qno, 1.0, Indices{index::e, index::E}, S_ijk_mijk_ein, Indices{index::e}, F_me_slice);

                    Tensor<double, 3> R_ijk_buffer_b("R_ijk_buffer_b", n_qno_[mijk], n_qno_[mijk], n_qno_[mijk]);
                    einsum(0.0, Indices{index::a, index::b, index::c}, &R_ijk_buffer_b, 1.0 / 3, Indices{index::e, index::a, index::b, index::c}, 
                            alpha_mijk, Indices{index::e}, F_me_qno);
                    
                    R_ijk_buffer += matmul_3d_einsums(R_ijk_buffer_b, S_ijk_mijk, n_qno_[mijk], n_tno_[ijk]);

                    // Jiang and Matthews Eq. 6b 1/2 (ae|mf) alpha_{mijk}^{febc}
                    Tensor<double, 3> K_ledb_slice = K_ledb(m_ijk, All, All, All);
                    Tensor<double, 3> K_ledb_temp_a = matmul_3d_index(K_ledb_slice, S_ijk_mijk->transpose(), 0);
                    Tensor<double, 3> K_ledb_temp_b = matmul_3d_index(K_ledb_temp_a, S_ijk_mijk->transpose(), 1);

                    Tensor<double, 3> R_ijk_buffer_c("R_ijk_buffer_c", n_tno_[ijk], n_qno_[mijk], n_qno_[mijk]);
                    einsum(0.0, Indices{index::a, index::b, index::c}, &R_ijk_buffer_c, 1.0, Indices{index::f, index::e, index::a}, K_ledb_temp_b,
                            Indices{index::f, index::e, index::b, index::c}, alpha_mijk);
                    
                    Tensor<double, 3> R_ijk_buffer_d = matmul_3d_index(R_ijk_buffer_c, S_ijk_mijk, 1);
                    R_ijk_buffer += matmul_3d_index(R_ijk_buffer_d, S_ijk_mijk, 2);
                }

                for (int n_ijk = 0; n_ijk < nlmo_ijk; ++n_ijk) {
                    int n = lmotriplet_to_lmos_[ijk][n_ijk];
                    int mnjk_idx = m * std::pow(naocc, 3) + n * std::pow(naocc, 2) + j * naocc + k;
                    int mnjk = i_j_k_l_to_ijkl_.count(mnjk_idx) ? (i_j_k_l_to_ijkl_[mnjk_idx]) : -1;
                    if (mnjk == -1) continue;

                    // Jiang and Matthews Eq. 6c -1/2 (me|ni) alpha_{mnjk}^{eabc} -> O(N^{10})
                    auto S_ijk_mnjk = submatrix_rows_and_cols(*S_pao_, lmotriplet_to_paos_[ijk], lmoquadruplet_to_paos_[mnjk]);
                    S_ijk_mnjk = linalg::triplet(X_tno_[ijk], S_ijk_mnjk, X_qno_[mnjk], true, false, false);
                    Tensor<double, 2> S_ijk_mnjk_ein("S_ijk_mnjk_ein", n_tno_[ijk], n_qno_[mnjk]);
                    ::memcpy(S_ijk_mnjk_ein.data(), S_ijk_mnjk->get_pointer(), n_tno_[ijk] * n_qno_[mnjk] * sizeof(double));

                    Tensor<double, 4> T_mnjk = quadruples_permuter(T_iajbkcld_[mnjk], m, n, j, k);
                    Tensor<double, 4> alpha_mnjk = alpha_ijkl_helper(T_mnjk);

                    Tensor<double, 1> K_meni_slice = K_menj_list[i_idx](All, m_ijk, n_ijk);
                    Tensor<double, 1> K_meni_qno("K_meni_qno", n_qno_[mnjk]);
                    einsum(0.0, Indices{index::E}, &K_meni_qno, 1.0, Indices{index::e, index::E}, S_ijk_mnjk_ein, Indices{index::e}, K_meni_slice);

                    Tensor<double, 3> R_ijk_buffer_e("R_ijk_buffer_e", n_qno_[mnjk], n_qno_[mnjk], n_qno_[mnjk]);
                    einsum(0.0, Indices{index::a, index::b, index::c}, &R_ijk_buffer_e, -1.0, Indices{index::e, index::a, index::b, index::c}, 
                            alpha_mnjk, Indices{index::e}, K_meni_qno);

                    R_ijk_buffer += matmul_3d_einsums(R_ijk_buffer_e, S_ijk_mnjk, n_qno_[mnjk], n_tno_[ijk]);
                } // end n_ijk
            } // end m_ijk
            R_ijk_cont += triples_permuter_einsums(R_ijk_buffer, i_idx, j_idx, k_idx); // Lie algebra rules for short perms allow for this
        } // end perm_idx

        C_DAXPY(ntno_ijk * ntno_ijk * ntno_ijk, 1.0, R_ijk_cont.data(), 1, R_iajbkc[ijk]->get_pointer(), 1);

        if (disk_ints_) {
            q_ov_[ijk] = Tensor<double, 3>(q_ov_[ijk].name(), 0, 0, 0);
            q_vv_[ijk] = Tensor<double, 3>(q_vv_[ijk].name(), 0, 0, 0);
        }
    }
}

// Jiang and Matthews Eq. 28b
std::vector<Tensor<double, 4>> DLPNOCCSDTQ::alpha_L_contribution() {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // 0.5 (me|nf) alpha_{nijk}^{fabe}
    std::vector<Tensor<double, 4>> L_alpha_list(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int naux_ij = lmopair_to_ribfs_[ij].size();
        int npno_ij = n_pno_[ij];

        L_alpha_list[ij] = Tensor<double, 4>("L_alpha", nlmo_ij, nlmo_ij, npno_ij, npno_ij); // (k, m, a, b)
        L_alpha_list[ij].zero();

        std::vector<SharedMatrix> q_ov_ij = QIA_PNO(ij);

        Tensor<double, 4> g_menf_t("g_menf_t", nlmo_ij, nlmo_ij, npno_ij, npno_ij);
        g_menf_t.zero();

        for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
            int m = lmopair_to_lmos_[ij][m_ij];
            for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                int n = lmopair_to_lmos_[ij][n_ij];
                for (int e_ij = 0; e_ij < npno_ij; ++e_ij) {
                    for (int f_ij = 0; f_ij < npno_ij; ++f_ij) {
                        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                            g_menf_t(m_ij, n_ij, e_ij, f_ij) += (*q_ov_ij[q_ij])(m_ij, e_ij) * (*q_ov_ij[q_ij])(n_ij, f_ij);
                        } // end q_ij
                    } // end f_ij
                } // end e_ij
            } // end n_ij
        } // end m_ij
        
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];

            Tensor<double, 3> L_alpha_temp("L_alpha_temp", nlmo_ij, npno_ij, npno_ij);
            L_alpha_temp.zero();

            for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                int n = lmopair_to_lmos_[ij][n_ij];

                int nijk_idx = n * std::pow(naocc, 3) + i * std::pow(naocc, 2) + j * naocc + k;
                int nijk = i_j_k_l_to_ijkl_.count(nijk_idx) ? i_j_k_l_to_ijkl_[nijk_idx] : -1;

                if (nijk == -1) continue;

                auto S_nijk_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[nijk], lmopair_to_paos_[ij]);
                S_nijk_ij = linalg::triplet(X_qno_[nijk], S_nijk_ij, X_pno_[ij], true, false, false);

                Tensor<double, 4> T_nijk = matmul_4d(quadruples_permuter(T_iajbkcld_[nijk], n, i, j, k), S_nijk_ij->transpose(), n_qno_[nijk], n_pno_[ij]);
                Tensor<double, 4> alpha_nijk = alpha_ijkl_helper(T_nijk);
                Tensor<double, 3> g_menf_t_slice = g_menf_t(n_ij, All, All, All); // (m, f, e)

                permute(Indices{index::f, index::e, index::a, index::b}, &T_nijk, Indices{index::f, index::a, index::b, index::e}, alpha_nijk);
                
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_alpha_temp, 0.5, Indices{index::m, index::f, index::e}, g_menf_t_slice,
                        Indices{index::f, index::e, index::a, index::b}, T_nijk);

            } // end n_ij

            L_alpha_list[ij](k_ij, All, All, All) = L_alpha_temp;

        } // end k_ij
    }

    return L_alpha_list;

}

// Jiang and Matthews Eq. 29b
std::vector<Tensor<double, 4>> DLPNOCCSDTQ::alpha_M_contribution() {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // -0.5 (me|nf) alpha_{nmjk}^{fabc}
    std::vector<Tensor<double, 4>> M_alpha_list(n_lmo_pairs);

#pragma omp parallel for
    for (int jk = 0; jk < n_lmo_pairs; ++jk) {
        auto &[j, k] = ij_to_i_j_[jk];
        int nlmo_jk = lmopair_to_lmos_[jk].size();
        int naux_jk = lmopair_to_ribfs_[jk].size();
        int npno_jk = n_pno_[jk];

        M_alpha_list[jk] = Tensor<double, 4>("M_alpha", npno_jk, npno_jk, npno_jk, npno_jk);
        M_alpha_list[jk].zero();

        std::vector<SharedMatrix> q_ov_jk = QIA_PNO(jk);

        Tensor<double, 4> g_menf_t("g_menf_t", nlmo_jk, nlmo_jk, npno_jk, npno_jk);
        g_menf_t.zero();

        for (int n_jk = 0; n_jk < nlmo_jk; ++n_jk) {
            int n = lmopair_to_lmos_[jk][n_jk];
            for (int m_jk = 0; m_jk < nlmo_jk; ++m_jk) {
                int m = lmopair_to_lmos_[jk][m_jk];
                for (int e_jk = 0; e_jk < n_pno_[jk]; ++e_jk) {
                    for (int f_jk = 0; f_jk < n_pno_[jk]; ++f_jk) {
                        for (int q_jk = 0; q_jk < naux_jk; ++q_jk) {
                            g_menf_t(n_jk, m_jk, f_jk, e_jk) += (*q_ov_jk[q_jk])(n_jk, f_jk) * (*q_ov_jk[q_jk])(m_jk, e_jk);
                        } // end q_jk
                    } // end f_jk
                } // end e_jk
            } // end m_jk
        } // end n_jk

        for (int n_jk = 0; n_jk < nlmo_jk; ++n_jk) {
            int n = lmopair_to_lmos_[jk][n_jk];
            for (int m_jk = 0; m_jk < nlmo_jk; ++m_jk) {
                int m = lmopair_to_lmos_[jk][m_jk];
                int nmjk_idx = n * std::pow(naocc, 3) + m * std::pow(naocc, 2) + j * naocc + k;
                int nmjk = i_j_k_l_to_ijkl_.count(nmjk_idx) ? i_j_k_l_to_ijkl_[nmjk_idx] : -1;

                if (nmjk == -1) continue;

                auto S_nmjk_jk = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[nmjk], lmopair_to_paos_[jk]);
                S_nmjk_jk = linalg::triplet(X_qno_[nmjk], S_nmjk_jk, X_pno_[jk], true, false, false);

                Tensor<double, 4> T_nmjk = matmul_4d(quadruples_permuter(T_iajbkcld_[nmjk], n, m, j, k), S_nmjk_jk->transpose(), n_qno_[nmjk], n_pno_[jk]);
                Tensor<double, 4> alpha_nmjk = alpha_ijkl_helper(T_nmjk);
                Tensor<double, 2> g_menf_t_slice = g_menf_t(n_jk, m_jk, All, All);

                einsum(1.0, Indices{index::e, index::a, index::b, index::c}, &M_alpha_list[jk], -0.5, Indices{index::f, index::e}, g_menf_t_slice,
                        Indices{index::f, index::a, index::b, index::c}, alpha_nmjk);
            } // end n_jk
        } // end m_jk
    } // end jk

    return M_alpha_list;
}

void DLPNOCCSDTQ::form_T_mnkl() {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // Stored as kl, mn
    T_mnkl_list_.resize(n_lmo_pairs);

    // Loop over all pairs
#pragma omp parallel for schedule(dynamic, 1)
    for (int kl = 0; kl < n_lmo_pairs; ++kl) {
        auto &[k, l] = ij_to_i_j_[kl];
        if (k > l) continue;

        // number of LMOs in the quadruplet domain
        const int nlmo_kl = lmopair_to_lmos_[kl].size();
        // number of PNOs in the quadruplet domain
        const int npno_kl = n_pno_[kl];

        T_mnkl_list_[kl].resize(n_lmo_pairs);

        for (int m_kl = 0; m_kl < nlmo_kl; ++m_kl) {
            int m = lmopair_to_lmos_[kl][m_kl];
            for (int n_kl = 0; n_kl < nlmo_kl; ++n_kl) {
                int n = lmopair_to_lmos_[kl][n_kl];
                int mn = i_j_to_ij_[m][n];
                if (m > n) continue;

                int mnkl_idx = m * std::pow(naocc, 3) + n * std::pow(naocc, 2) + k * naocc + l;
                int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;

                if (mn == -1 || mnkl == -1) continue;
                
                auto S_mnkl_kl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[mnkl], lmopair_to_paos_ext_[kl]);
                S_mnkl_kl = linalg::triplet(X_qno_[mnkl], S_mnkl_kl, X_pno_ext_[kl], true, false, false);

                T_mnkl_list_[kl][mn] = matmul_4d(quadruples_permuter(T_iajbkcld_[mnkl], m, n, k, l), 
                        S_mnkl_kl->transpose(), n_qno_[mnkl], n_pno_ext_[kl]);
            } // end n_kl
        } // end m_kl
    } // end kl
}

void DLPNOCCSDTQ::compute_R_iajbkcld(std::vector<Tensor<double, 4>>& R_iajbkcld) {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // Index orders corresponding to the 24 occupied-orbital permutations.
    auto einsum_indices = std::make_tuple(Indices{a, b, c, d}, Indices{a, b, d, c}, Indices{a, c, b, d}, Indices{a, c, d, b}, 
        Indices{a, d, b, c}, Indices{a, d, c, b}, Indices{b, a, c, d}, Indices{b, a, d, c}, Indices{b, c, a, d}, Indices{b, c, d, a}, 
        Indices{b, d, a, c}, Indices{b, d, c, a}, Indices{c, a, b, d}, Indices{c, a, d, b}, Indices{c, b, a, d}, Indices{c, b, d, a}, 
        Indices{c, d, a, b}, Indices{c, d, b, a}, Indices{d, a, b, c}, Indices{d, a, c, b}, Indices{d, b, a, c}, Indices{d, b, c, a}, 
        Indices{d, c, a, b}, Indices{d, c, b, a});

    // Compute expensive alpha L contribution
    std::vector<Tensor<double, 4>> L_alpha_list = alpha_L_contribution();
    // Compute expensive alpha M contribution
    std::vector<Tensor<double, 4>> M_alpha_list = alpha_M_contribution();

// Loop over unique quadruplets
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
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

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // Bookkeeping for permutations of ijk (ijkl choose ijk)
        std::unordered_map<int, int> ijk_to_ijkl_perm_idx;

        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            int ijk_idx = i * std::pow(naocc, 2) + j * naocc + k;

            ijk_to_ijkl_perm_idx[ijk_idx] = perm_idx;
        });

        // This variable name sounds mean :(
        std::array<std::tuple<int, int, int>, 4> exclusion_list = {std::make_tuple(j, k, l), 
            std::make_tuple(i, k, l), std::make_tuple(j, i, l), std::make_tuple(j, k, i)};

        std::unordered_set<int> computed_perms;

        // (T1-DRESS) NECESSARY FOCK MATRIX ELEMENTS AND INTEGRALS

        // List of singles amplitudes projected onto ijkl space
        std::array<Tensor<double, 1>, 4> T_i_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            int i_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), i) - lmoquadruplet_to_lmos_[ijkl].begin();
            T_i_list[i_idx] = Tensor<double, 1>("T_i", nqno_ijkl);
            ::memcpy(T_i_list[i_idx].data(), &(T_n_ijkl_[ijkl])(i_ijkl, 0), nqno_ijkl * sizeof(double));
        } // end i_idx

        // Read in expensive integrals
        // (Q | m a) integrals
        Tensor<double, 3> q_ov = QIA_QNO(ijkl);
        // (Q | a b) integrals
        Tensor<double, 3> q_vv = QAB_QNO(ijkl);

        // T1-dress q_io and q_iv (check represents a T1-dressed index)
        std::array<Tensor<double, 2>, 4> q_io_t1_list; // (Q | m \check{i}) = (Q | m i) + (Q | m a) t_{i}^{a}
        std::array<Tensor<double, 2>, 4> q_iv_t1_list; // (Q | \check{a} \check{i}) = (Q | a i) - t_{m}^{a} (Q | m i) + (Q | a b) t_{i}^{b} - t_{m}^{a} (Q | m b) t_{i}^{b}

        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            q_io_t1_list[i_idx] = q_io_list_[ijkl][i_idx]; // (Q | m i)
            einsum(1.0, Indices{index::Q, index::m}, &q_io_t1_list[i_idx], 1.0, Indices{index::Q, index::m, index::a}, q_ov, Indices{index::a}, T_i_list[i_idx]);

            q_iv_t1_list[i_idx] = q_iv_list_[ijkl][i_idx]; // (Q | a i)
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], -1.0, Indices{index::Q, index::m}, q_io_list_[ijkl][i_idx], Indices{index::m, index::a}, T_n_ijkl_[ijkl]);
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], 1.0, Indices{index::Q, index::a, index::b}, q_vv, Indices{index::b}, T_i_list[i_idx]);
            Tensor<double, 2> q_iv_t1_temp("q_iv_t1_temp", naux_ijkl, nlmo_ijkl);
            einsum(0.0, Indices{index::Q, index::m}, &q_iv_t1_temp, 1.0, Indices{index::Q, index::m, index::b}, q_ov, Indices{index::b}, T_i_list[i_idx]);
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], -1.0, Indices{index::Q, index::m}, q_iv_t1_temp, Indices{index::m, index::a}, T_n_ijkl_[ijkl]);
        }

        // Store q_ov with transposed orbital indices so subsequent contractions can use GEMM kernels.
        Tensor<double, 3> q_vo("q_vo", naux_ijkl, nqno_ijkl, nlmo_ijkl);
        permute(Indices{index::Q, index::a, index::m}, &q_vo, Indices{index::Q, index::m, index::a}, q_ov);

        // T1-dress the second virtual index, opposite to the paper's convention, to reduce contraction cost.
        Tensor<double, 3> q_vv_t1 = q_vv; // (Q | a \check{b}) = (Q | a b) - (Q | a m) t_{m}^{b} 
        einsum(1.0, Indices{index::Q, index::a, index::b}, &q_vv_t1, -1.0, Indices{index::Q, index::a, index::m}, q_vo, Indices{index::m, index::b}, T_n_ijkl_[ijkl]);

        // Build the T1-dressed Fock intermediates.

        // Gamma_Q is used universally for J-like contractions
        Tensor<double, 1> gamma_Q("gamma_Q", naux_ijkl);
        einsum(0.0, Indices{index::Q}, &gamma_Q, 1.0, Indices{index::Q, index::m, index::e}, q_ov, Indices{index::m, index::e}, T_n_ijkl_[ijkl]);

        // F_me (this is scoped to ensure that the intermediate tensors are not persistent in memory)
        Tensor<double, 2> F_me("F_me", nlmo_ijkl, nqno_ijkl); {
            // J contractions
            einsum(0.0, Indices{index::m, index::e}, &F_me, 2.0, Indices{index::Q, index::m, index::e}, q_ov, Indices{index::Q}, gamma_Q);
            
            // K contractions (rc|ks)t_{k}^{c} ... (mf|ne) t_{n}^{f}
            Tensor<double, 3> F_me_K_temp("F_me_K_temp", naux_ijkl, nlmo_ijkl, nlmo_ijkl);
            einsum(0.0, Indices{index::Q, index::m, index::n}, &F_me_K_temp, 1.0, Indices{index::Q, index::m, index::f}, q_ov, Indices{index::n, index::f}, T_n_ijkl_[ijkl]);
            Tensor<double, 3> F_me_K_temp2("F_me_K_temp2", naux_ijkl, nlmo_ijkl, nlmo_ijkl);
            permute(Indices{index::Q, index::n, index::m}, &F_me_K_temp2, Indices{index::Q, index::m, index::n}, F_me_K_temp);
            einsum(1.0, Indices{index::m, index::e}, &F_me, -1.0, Indices{index::Q, index::n, index::m}, F_me_K_temp2, Indices{index::Q, index::n, index::e}, q_ov);
        }

        // F_mi (over i, j, k, l)
        std::array<Tensor<double, 1>, 4> F_mi_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];

            // F_lmo (non-T1 contribution)
            F_mi_list[i_idx] = Tensor<double, 1>("F_mi", nlmo_ijkl);
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                F_mi_list[i_idx](m_ijkl) = (*F_lmo_)(m, i);
            }

            einsum(1.0, Indices{index::m}, &F_mi_list[i_idx], 2.0, Indices{index::Q, index::m}, q_io_list_[ijkl][i_idx], Indices{index::Q}, gamma_Q);

            // K contractions (rc|ks)t_{k}^{c} ... (mf|ni) t_{n}^{f}
            Tensor<double, 2> F_mi_K_temp("F_mi_K_temp", naux_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::Q, index::f}, &F_mi_K_temp, 1.0, Indices{index::Q, index::n}, q_io_list_[ijkl][i_idx], Indices{index::n, index::f}, T_n_ijkl_[ijkl]);
            einsum(1.0, Indices{index::m}, &F_mi_list[i_idx], -1.0, Indices{index::Q, index::f, index::m}, q_vo, Indices{index::Q, index::f}, F_mi_K_temp);

            // Add F_me contribution
            einsum(1.0, Indices{index::m}, &F_mi_list[i_idx], 1.0, Indices{index::m, index::e}, F_me, Indices{index::e}, T_i_list[i_idx]);
        }

        // => F_ae <= //
        Tensor<double, 2> F_ae("F_ae", nqno_ijkl, nqno_ijkl); {
            F_ae.zero();
            // e_qno (non-t1 contribution)
            for (int a_ijkl = 0; a_ijkl < nqno_ijkl; ++a_ijkl) {
                F_ae(a_ijkl, a_ijkl) = (*e_qno_[ijkl])(a_ijkl);
            }

            // J contribution
            einsum(1.0, Indices{index::a, index::e}, &F_ae, 2.0, Indices{index::Q, index::a, index::e}, q_vv, Indices{index::Q}, gamma_Q);
            // K contribution (rc|ks)t_{k}^{c} ... (af|ne) t_{n}^{f}
            Tensor<double, 3> F_ae_K_temp("F_ae_K_temp", naux_ijkl, nqno_ijkl, nlmo_ijkl);
            einsum(0.0, Indices{index::Q, index::a, index::n}, &F_ae_K_temp, 1.0, Indices{index::Q, index::a, index::f}, q_vv, Indices{index::n, index::f}, T_n_ijkl_[ijkl]);
            Tensor<double, 3> F_ae_K_temp2("F_ae_K_temp2", naux_ijkl, nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::Q, index::n, index::a}, &F_ae_K_temp2, Indices{index::Q, index::a, index::n}, F_ae_K_temp);
            einsum(1.0, Indices{index::a, index::e}, &F_ae, -1.0, Indices{index::Q, index::n, index::a}, F_ae_K_temp2, Indices{index::Q, index::n, index::e}, q_ov);

            // Add the F_me contribution to F_ae
            einsum(1.0, Indices{index::a, index::e}, &F_ae, -1.0, Indices{index::m, index::a}, T_n_ijkl_[ijkl], Indices{index::m, index::e}, F_me);
        }

        // Amplitude intermediates
        std::array<Tensor<double, 2>, 16> T_ij_list; // (Project T_ij amplitudes from PNO space of ij to QNO space of ijkl)
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j];

                T_ij_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("T_ij", nqno_ijkl, nqno_ijkl);
                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);

                auto T_ij = linalg::triplet(S_ijkl_ij, T_iajb_[ij], S_ijkl_ij, false, false, true);
                ::memcpy(T_ij_list[i_idx * FOUR + j_idx].data(), T_ij->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            } // end j_idx
        } // end i_idx

        std::array<Tensor<double, 3>, 4> T_mi_list; // Project T_mi amplitudes from PNO space of mi to QNO space of ijkl
        std::array<Tensor<double, 3>, 4> U_mi_list; // Project U_mi amplitudes from PNO space of mi to QNO space of ijkl
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            T_mi_list[i_idx] = Tensor<double, 3>("T_mi", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            U_mi_list[i_idx] = Tensor<double, 3>("U_mi", nlmo_ijkl, nqno_ijkl, nqno_ijkl);

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                int mi = i_j_to_ij_[m][i];
                auto S_ijkl_mi = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[mi]);
                S_ijkl_mi = linalg::triplet(X_qno_[ijkl], S_ijkl_mi, X_pno_[mi], true, false, false);

                auto T_mi = linalg::triplet(S_ijkl_mi, T_iajb_[mi], S_ijkl_mi, false, false, true);
                auto U_mi = linalg::triplet(S_ijkl_mi, Tt_iajb_[mi], S_ijkl_mi, false, false, true);

                ::memcpy(&T_mi_list[i_idx](m_ijkl, 0, 0), T_mi->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
                ::memcpy(&U_mi_list[i_idx](m_ijkl, 0, 0), U_mi->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
            } // end m_ijkl
        } // end i_idx

        Tensor<double, 4> T_mn("T_mn", nlmo_ijkl, nlmo_ijkl, nqno_ijkl, nqno_ijkl); {
            // Project T_mn amplitudes from PNO space of mn to QNO space of ijkl
            T_mn.zero();
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                    int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];
                    int mn = i_j_to_ij_[m][n];
                    if (mn == -1) continue;

                    auto S_ijkl_mn = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[mn]);
                    S_ijkl_mn = linalg::triplet(X_qno_[ijkl], S_ijkl_mn, X_pno_[mn], true, false, false);

                    auto T_mn_ijkl = linalg::triplet(S_ijkl_mn, T_iajb_[mn], S_ijkl_mn, false, false, true);
                    ::memcpy(&T_mn(m_ijkl, n_ijkl, 0, 0), T_mn_ijkl->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));
                }
            }
        }

        std::array<Tensor<double, 5>, 4> T_mni_list; // Project T_mni amplitudes from TNO space of mni to QNO space of ijkl
        // This contraction has O(N^{10}) worst-case cost.
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            T_mni_list[i_idx] = Tensor<double, 5>("T_mni", nlmo_ijkl, nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            T_mni_list[i_idx].zero();

            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                    int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];
                    int mni_idx = m * std::pow(naocc, 2) + n * naocc + i;
                    int mni = (i_j_k_to_ijk_.count(mni_idx)) ? i_j_k_to_ijk_[mni_idx] : -1;
                    if (mni == -1) continue;

                    auto S_ijkl_mni = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[mni]);
                    S_ijkl_mni = linalg::triplet(X_qno_[ijkl], S_ijkl_mni, X_tno_[mni], true, false, false);
                    auto T_mni = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[mni], m, n, i), 
                                                    S_ijkl_mni, n_tno_[mni], n_qno_[ijkl]); // O(N^{10}) // (a, b, c) (c, C)
                    
                    ::memcpy(&T_mni_list[i_idx](m_ijkl, n_ijkl, 0, 0, 0), T_mni.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                } // end n_ijkl
            } // end m_ijkl
        } // end i_idx

        // Project T_mij amplitudes from TNO space of mij to QNO space of ijkl
        std::array<Tensor<double, 4>, 16> T_mij_list;
        std::array<Tensor<double, 4>, 16> Z_mij_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                int j = ijkl_list[j_idx];
                T_mij_list[i_idx * FOUR + j_idx] = Tensor<double, 4>("T_mij", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                Z_mij_list[i_idx * FOUR + j_idx] = Tensor<double, 4>("Z_mij", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                T_mij_list[i_idx * FOUR + j_idx].zero();
                Z_mij_list[i_idx * FOUR + j_idx].zero();

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    int mij_idx = m * std::pow(naocc, 2) + i * naocc + j;
                    int mij = i_j_k_to_ijk_.count(mij_idx) ? i_j_k_to_ijk_[mij_idx] : -1;
                    if (mij == -1) continue;

                    auto S_ijkl_mij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[mij]);
                    S_ijkl_mij = linalg::triplet(X_qno_[ijkl], S_ijkl_mij, X_tno_[mij], true, false, false);

                    Tensor<double, 3> T_mij = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[mij], m, i, j), 
                                                S_ijkl_mij, n_tno_[mij], n_qno_[ijkl]);
                    // Z_mij^{abc} = 2 T_mij^{abc} - T_mij^{bac} - T_mij^{cba}.
                    Tensor<double, 3> T_mij_bac("T_mij_bac", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::b, index::a, index::c}, &T_mij_bac, Indices{index::a, index::b, index::c}, T_mij);
                    Tensor<double, 3> T_mij_cba("T_mij_cba", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::c, index::b, index::a}, &T_mij_cba, Indices{index::a, index::b, index::c}, T_mij);

                    Tensor<double, 3> Z_mij = T_mij;
                    Z_mij *= 2.0;
                    Z_mij -= T_mij_bac;
                    Z_mij -= T_mij_cba;

                    ::memcpy(&T_mij_list[i_idx * FOUR + j_idx](m_ijkl, 0, 0, 0), T_mij.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                    ::memcpy(&Z_mij_list[i_idx * FOUR + j_idx](m_ijkl, 0, 0, 0), Z_mij.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                }
            } // end j_idx
        } // end i_idx

        // Quadruples amplitude intermediates
        // (i -> jkl), (j -> ikl), (k -> ijl), (l -> ijk)
        std::array<Tensor<double, 5>, 4> T_nijk_exclusion_list;
        std::array<Tensor<double, 5>, 4> alpha_nijk_exclusion_list;
        std::array<Tensor<double, 5>, 4> T_nijk_exclusion_list_unsorted;

        for (int idx = 0; idx < FOUR; ++idx) {
            auto &[i, j, k] = exclusion_list[idx];
            
            T_nijk_exclusion_list[idx] = Tensor<double, 5>("T_nijk", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            T_nijk_exclusion_list[idx].zero();

            alpha_nijk_exclusion_list[idx] = Tensor<double, 5>("alpha_nijk", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            alpha_nijk_exclusion_list[idx].zero();

            T_nijk_exclusion_list_unsorted[idx] = Tensor<double, 5>("T_nijk_unsorted", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            T_nijk_exclusion_list_unsorted[idx].zero();

            for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];
                int nijk_idx = n * std::pow(naocc, 3) + i * std::pow(naocc, 2) + j * naocc + k;
                int nijk = i_j_k_l_to_ijkl_.count(nijk_idx) ? (i_j_k_l_to_ijkl_[nijk_idx]) : -1;
                if (nijk == -1) continue;

                auto S_ijkl_nijk = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmoquadruplet_to_paos_[nijk]);
                S_ijkl_nijk = linalg::triplet(X_qno_[ijkl], S_ijkl_nijk, X_qno_[nijk], true, false, false);

                Tensor<double, 4> T_nijk_unsorted = matmul_4d(T_iajbkcld_[nijk], S_ijkl_nijk, n_qno_[nijk], n_qno_[ijkl]);
                ::memcpy(&T_nijk_exclusion_list_unsorted[idx](n_ijkl, 0, 0, 0, 0), T_nijk_unsorted.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                
                Tensor<double, 4> T_nijk = quadruples_permuter(T_nijk_unsorted, n, i, j, k);
                ::memcpy(&T_nijk_exclusion_list[idx](n_ijkl, 0, 0, 0, 0), T_nijk.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));

                Tensor<double, 4> alpha_nijk = alpha_ijkl_helper(T_nijk);
                ::memcpy(&alpha_nijk_exclusion_list[idx](n_ijkl, 0, 0, 0, 0), alpha_nijk.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                
            } // end n_ijkl
        } // end i_idx

        // Quadruples amplitude intermediates

        // => Form integrals <= //

        // (ov|oo) integrals
        // B^{Q}_{me} B^{Q}_{nj}
        std::array<Tensor<double, 3>, 4> g_menj_list; // Stored as (e, m, n)
        std::array<Tensor<double, 3>, 4> g_menj_t_list; // Stored as (m, e, n)

        // 2 B^{Q}_{me} B^{Q}_{nj} - B^{Q}_{ne} B^{Q}_{mj}
        std::array<Tensor<double, 3>, 4> L_menj_list; // Stored as (e, m, n)
        std::array<Tensor<double, 3>, 4> L_menj_t_list; // Stored as (m, e, n)

        for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
            g_menj_list[j_idx] = Tensor<double, 3>("(me|nj)", nqno_ijkl, nlmo_ijkl, nlmo_ijkl); // (e, m, n)
            einsum(0.0, Indices{index::e, index::m, index::n}, &g_menj_list[j_idx], 1.0, Indices{index::Q, index::e, index::m}, q_vo,
                    Indices{index::Q, index::n}, q_io_t1_list[j_idx]);
            g_menj_t_list[j_idx] = Tensor<double, 3>("g_menj_t", nlmo_ijkl, nqno_ijkl, nlmo_ijkl);
            permute(Indices{index::m, index::e, index::n}, &g_menj_t_list[j_idx], Indices{index::e, index::m, index::n}, g_menj_list[j_idx]);

            Tensor<double, 3> L_menj_temp("L_menj", nqno_ijkl, nlmo_ijkl, nlmo_ijkl);
            permute(Indices{index::e, index::n, index::m}, &L_menj_temp, Indices{index::e, index::m, index::n}, g_menj_list[j_idx]);

            L_menj_list[j_idx] = g_menj_list[j_idx];
            L_menj_list[j_idx] *= 2.0;
            L_menj_list[j_idx] -= L_menj_temp;

            L_menj_t_list[j_idx] = Tensor<double, 3>("L_menj_t", nlmo_ijkl, nqno_ijkl, nlmo_ijkl);
            permute(Indices{index::m, index::e, index::n}, &L_menj_t_list[j_idx], Indices{index::e, index::m, index::n}, L_menj_list[j_idx]);
        }

        // (ov|ov) integrals
        // B^{Q}_{me} B^{Q}_{nf}
        Tensor<double, 4> g_menf("(me|nf)", nlmo_ijkl, nqno_ijkl, nlmo_ijkl, nqno_ijkl); // Stored as (m, e, n, f)
        einsum(0.0, Indices{index::m, index::e, index::n, index::f}, &g_menf, 1.0, Indices{index::Q, index::m, index::e}, q_ov,
                Indices{index::Q, index::n, index::f}, q_ov);

        // 2 B^{Q}_{me} B^{Q}_{nf} - B^{Q}_{mf} B^{Q}_{ne}
        Tensor<double, 4> L_menf = g_menf; {
            Tensor<double, 4> L_temp("L_temp", nlmo_ijkl, nqno_ijkl, nlmo_ijkl, nqno_ijkl); // Stored as (m, e, n, f)
            permute(Indices{index::m, index::f, index::n, index::e}, &L_temp, Indices{index::m, index::e, index::n, index::f}, g_menf);
            L_menf *= 2.0;
            L_menf -= L_temp;
        }
        
        Tensor<double, 4> g_menf_t("<mn|ef>", nlmo_ijkl, nlmo_ijkl, nqno_ijkl, nqno_ijkl); // Stored as (m, n, e, f)... Physicist's notation
        permute(Indices{index::m, index::n, index::e, index::f}, &g_menf_t, Indices{index::m, index::e, index::n, index::f}, g_menf);
        Tensor<double, 4> L_menf_t("L_mnef_t", nlmo_ijkl, nlmo_ijkl, nqno_ijkl, nqno_ijkl); // Stored as (m, n, e, f)
        permute(Indices{index::m, index::n, index::e, index::f}, &L_menf_t, Indices{index::m, index::e, index::n, index::f}, L_menf);

        // (ov|vo) integrals

        // B^{Q}_{me} B^{Q}_{ai}
        std::array<Tensor<double, 3>, 4> g_meai_list; // Stored as (m, e, a)
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            g_meai_list[i_idx] = Tensor<double, 3>("(me|ai)", nlmo_ijkl, nqno_ijkl, nqno_ijkl); // (m, e, a)
            einsum(0.0, Indices{index::m, index::e, index::a}, &g_meai_list[i_idx], 1.0, Indices{index::Q, index::m, index::e},
                    q_ov, Indices{index::Q, index::a}, q_iv_t1_list[i_idx]);
        }

        // (oo|vv) integrals

        // B^{Q}_{mi} B^{Q}_{ae}
        std::array<Tensor<double, 3>, 4> g_miae_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            g_miae_list[i_idx] = Tensor<double, 3>("(mi|ae)", nlmo_ijkl, nqno_ijkl, nqno_ijkl); // (m, e, a)
            einsum(0.0, Indices{index::m, index::e, index::a}, &g_miae_list[i_idx], 1.0, Indices{index::Q, index::m}, q_io_t1_list[i_idx],
                    Indices{index::Q, index::e, index::a}, q_vv_t1);
        }

        // (ov|vv) like terms (different transposes, t and u are generated to take advantage of GEMM)

        // B^{Q}_{mf} B^{Q}_{ae}
        Tensor<double, 4> g_mfae = Tensor<double, 4>("(mf|ae)", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
        einsum(0.0, Indices{index::m, index::f, index::e, index::a}, &g_mfae, 1.0, Indices{index::Q, index::m, index::f}, q_ov,
                Indices{index::Q, index::e, index::a}, q_vv_t1); // Stored as (m, f, e, a)
        Tensor<double, 4> g_mfae_t("g_mfae_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl); // Stored as (m, e, f, a)
        permute(Indices{index::m, index::e, index::f, index::a}, &g_mfae_t, Indices{index::m, index::f, index::e, index::a}, g_mfae);
        Tensor<double, 4> g_mfae_u("g_mfae_u", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl); // Stored as (m, a, f, e)
        permute(Indices{index::m, index::a, index::f, index::e}, &g_mfae_u, Indices{index::m, index::f, index::e, index::a}, g_mfae);

        // 2 B^{Q}_{mf} B^{Q}_{ae} - B^{Q}_{me} B^{Q}_{af}
        Tensor<double, 4> L_mfae = g_mfae; // (m, f, e, a)
        L_mfae *= 2.0;
        L_mfae -= g_mfae_t; // (m, e, f, a)
        Tensor<double, 4> L_mfae_t = Tensor<double, 4>("L_mfae_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl); // (m, f, a, e)
        permute(Indices{index::m, index::f, index::a, index::e}, &L_mfae_t, Indices{index::m, index::f, index::e, index::a}, L_mfae);

        // A_{ej}^{ab} (dimensions: 4 * (e, a, b)) (Jiang and Matthews Eq. 7)
        std::array<Tensor<double, 3>, 4> A_ejab_list;
        for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
            int j = ijkl_list[j_idx];

            // Jiang Eq. 7a [(ae|bj)]
            A_ejab_list[j_idx] = Tensor<double, 3>("A_ejab", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], 1.0, Indices{index::Q, index::e, index::a}, 
                    q_vv_t1, Indices{index::Q, index::b}, q_iv_t1_list[j_idx]);

            // Jiang Eq. 7b [(me|nj)T_{mn}^{ab}]
            einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], 1.0, Indices{index::e, index::m, index::n}, g_menj_list[j_idx],
                    Indices{index::m, index::n, index::a, index::b}, T_mn);

            // Jiang Eq. 7c 1/2 * [2(mf|ae) - (me|af)]U_{mj}^{fb}
            einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], 0.5, Indices{index::m, index::f, index::e, index::a},
                    L_mfae, Indices{index::m, index::f, index::b}, U_mi_list[j_idx]);

            // Jiang Eq. 7d -(1/2 + P_{ab})(me|af)T_{mj}^{bf}
            // Form the direct and permuted contributions to Jiang Eq. 7d.
            Tensor<double, 3> A_eq7d_direct("A_eq7d_direct", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 3> A_eq7d_permuted("A_eq7d_permuted", nqno_ijkl, nqno_ijkl, nqno_ijkl);

            Tensor<double, 3> T_mi_t("T_mi_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::f, index::b}, &T_mi_t, Indices{index::m, index::b, index::f}, T_mi_list[j_idx]);

            einsum(0.0, Indices{index::e, index::a, index::b}, &A_eq7d_direct, 1.0, Indices{index::m, index::f, index::e, index::a},
                    g_mfae_t, Indices{index::m, index::f, index::b}, T_mi_t);
            permute(Indices{index::e, index::b, index::a}, &A_eq7d_permuted, Indices{index::e, index::a, index::b}, A_eq7d_direct);
            A_eq7d_direct *= 0.5;
            A_ejab_list[j_idx] -= A_eq7d_direct;
            A_ejab_list[j_idx] -= A_eq7d_permuted;

            // Jiang Eq. 7e -(me|nf) Z_{nmj}^{fab}
            // Z_{ijk}^{abc} = 2 T_{ijk}^{abc} - T_{ijk}^{bac} - T_{ijk}^{cba}
            // Form Z_{nmj} one triplet at a time to limit peak memory.

            for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    Tensor<double, 2> g_nm_slice = g_menf_t(n_ijkl, m_ijkl, All, All); // (f, e)

                    int nmj_idx = n * std::pow(naocc, 2) + m * naocc + j;
                    int nmj = (i_j_k_to_ijk_.count(nmj_idx)) ? i_j_k_to_ijk_[nmj_idx] : -1;

                    if (nmj != -1) {
                        auto S_ijkl_nmj = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[nmj]);
                        S_ijkl_nmj = linalg::triplet(X_qno_[ijkl], S_ijkl_nmj, X_tno_[nmj], true, false, false);
                        Tensor<double, 3> Z_nmj = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[nmj], n, m, j),
                                                                        S_ijkl_nmj, n_tno_[nmj], n_qno_[ijkl]);
                        // Z_nmj^{abc} = 2 T_nmj^{abc} - T_nmj^{bac} - T_nmj^{cba}.
                        Tensor<double, 3> T_nmj_bac("T_nmj_bac", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                        permute(Indices{index::b, index::a, index::c}, &T_nmj_bac, Indices{index::a, index::b, index::c}, Z_nmj);
                        Tensor<double, 3> T_nmj_cba("T_nmj_cba", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                        permute(Indices{index::c, index::b, index::a}, &T_nmj_cba, Indices{index::a, index::b, index::c}, Z_nmj);

                        Z_nmj *= 2.0;
                        Z_nmj -= T_nmj_bac;
                        Z_nmj -= T_nmj_cba;

                        // Accumulate the Eq. 7e contraction.
                        einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], -1.0, Indices{index::f, index::e},
                                g_nm_slice, Indices{index::f, index::a, index::b}, Z_nmj);
                    }
                } // end m_ijkl
            } // end n_ijkl

            // Jiang Eq. 7f -F_{me} T_{mj}^{ab}
            einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], -1.0, Indices{index::m, index::e}, F_me,
                    Indices{index::m, index::a, index::b}, T_mi_list[j_idx]);
        }

        // B_{ij}^{am} (dimensions: 16 * (m, a)) (Jiang and Matthews Eq. 9)
        std::array<Tensor<double, 2>, 16> B_ijam_list;

        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];

            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                int j = ijkl_list[j_idx];

                // Jiang Eq. 9a (ai|mj)
                B_ijam_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("B_ijam", nlmo_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::m}, q_io_t1_list[j_idx], 
                        Indices{index::Q, index::a}, q_iv_t1_list[i_idx]);

                // Jiang Eq. 9b (mf|ae)T_{ji}^{fe}
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::m, index::a, index::f, index::e},
                        g_mfae_u, Indices{index::f, index::e}, T_ij_list[j_idx * FOUR + i_idx]);

                // Jiang Eq. 9c 0.5 * [2(ne|mj) - (me|nj)]U_{ni}^{ea}
                Tensor<double, 3> U_mi_t("U_mi_t", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
                permute(Indices{index::e, index::n, index::a}, &U_mi_t, Indices{index::n, index::e, index::a}, U_mi_list[i_idx]);
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 0.5, Indices{index::e, index::n, index::m}, 
                        L_menj_list[j_idx], Indices{index::e, index::n, index::a}, U_mi_t);

                // Jiang Eq. 9d -(0.5 + P_{ij}) (me|nj)T_{ni}^{ae}
                Tensor<double, 3> T_mi_t("T_mi_t", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
                permute(Indices{index::a, index::n, index::e}, &T_mi_t, Indices{index::n, index::a, index::e}, T_mi_list[i_idx]);
                Tensor<double, 3> g_menj_u("g_menj_u", nlmo_ijkl, nlmo_ijkl, nqno_ijkl);
                permute(Indices{index::m, index::n, index::e}, &g_menj_u, Indices{index::e, index::m, index::n}, g_menj_list[j_idx]);
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], -0.5, Indices{index::m, index::n, index::e},
                        g_menj_u, Indices{index::a, index::n, index::e}, T_mi_t);

                Tensor<double, 3> T_mj_t("T_mj_t", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
                permute(Indices{index::a, index::n, index::e}, &T_mj_t, Indices{index::n, index::a, index::e}, T_mi_list[j_idx]);
                Tensor<double, 3> g_meni_u("g_meni_u", nlmo_ijkl, nlmo_ijkl, nqno_ijkl);
                permute(Indices{index::m, index::n, index::e}, &g_meni_u, Indices{index::e, index::m, index::n}, g_menj_list[i_idx]);
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], -1.0, Indices{index::m, index::n, index::e},
                        g_meni_u, Indices{index::a, index::n, index::e}, T_mj_t);

                // Jiang Eq. 9e F_{me} T_{ij}^{ae}
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::m, index::e},
                        F_me, Indices{index::a, index::e}, T_ij_list[i_idx * FOUR + j_idx]);
                        
                // Jiang Eq. 9f (me|nf) Z_{nij}^{fae}
                // Transpose Z_nij for the Eq. 9f contraction.
                Tensor<double, 4> Z_nij_transposed("Z_nij_transposed", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::n, index::e, index::f, index::a}, &Z_nij_transposed, Indices{index::n, index::f, index::a, index::e}, Z_mij_list[i_idx * FOUR + j_idx]);
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::m, index::n, index::e, index::f}, g_menf_t,
                        Indices{index::n, index::e, index::f, index::a}, Z_nij_transposed);
            } // end j_idx
        } // end i_idx

        // C_{ae} (dimensions: (a, e)) (Jiang and Matthews Eq. 11)
        Tensor<double, 2> C_ae = F_ae;
        einsum(1.0, Indices{index::a, index::e}, &C_ae, -1.0, Indices{index::n, index::m, index::f, index::a}, T_mn,
                Indices{index::n, index::m, index::f, index::e}, L_menf_t);

        // D_{mi} (dimensions: 4 * (m)) (Jiang and Matthews Equation 12)
        std::array<Tensor<double, 1>, 4> D_mi_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            D_mi_list[i_idx] = F_mi_list[i_idx];

            Tensor<double, 3> T_mi_t("T_mi_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::n, index::e, index::f}, &T_mi_t, Indices{index::n, index::f, index::e}, T_mi_list[i_idx]);
            einsum(1.0, Indices{index::m}, &D_mi_list[i_idx], 1.0, Indices{index::m, index::n, index::e, index::f}, L_menf_t, 
                    Indices{index::n, index::e, index::f}, T_mi_t);
        }

        // E_{ei}^{ma} (dimensions: 4 * (m, e, a)) (Jiang and Matthews Eq. 13)
        std::array<Tensor<double, 3>, 4> E_eima_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            // Eq. 13a 2(me|ai) - (mi|ae)
            E_eima_list[i_idx] = g_meai_list[i_idx];
            E_eima_list[i_idx] *= 2.0;
            E_eima_list[i_idx] -= g_miae_list[i_idx];

            // Eq. 13b [2(me|nf) - (mf|ne)]U_{ni}^{fa}
            einsum(1.0, Indices{index::m, index::e, index::a}, &E_eima_list[i_idx], 1.0, Indices{index::m, index::e, index::n, index::f}, L_menf,
                    Indices{index::n, index::f, index::a}, U_mi_list[i_idx]);
        }

        // F_{ie}^{ma} (dimensions: 4 * (m, e, a)) (Jiang and Matthews Eq. 15)
        std::array<Tensor<double, 3>, 4> F_iema_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            // Eq. 13a (mi|ae)
            F_iema_list[i_idx] = g_miae_list[i_idx];

            // Eq. 13b -[(mf|ne)]T_{ni}^{af} (transpose later)
            // Recover (mf|ne) from L_menf = 2.0 * (me|nf) - (mf|ne).
            Tensor<double, 3> T_mi_t("T_mi_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::n, index::f, index::a}, &T_mi_t, Indices{index::n, index::a, index::f}, T_mi_list[i_idx]);

            einsum(1.0, Indices{index::m, index::e, index::a}, &F_iema_list[i_idx], -2.0, Indices{index::m, index::e, index::n, index::f}, g_menf,
                    Indices{index::n, index::f, index::a}, T_mi_t);
            einsum(1.0, Indices{index::m, index::e, index::a}, &F_iema_list[i_idx], 1.0, Indices{index::m, index::e, index::n, index::f}, L_menf,
                    Indices{index::n, index::f, index::a}, T_mi_t);
        }

        // G_{ij}^{mn} (dimensions: 16 * (m, n)) (Jiang and Matthews Eq. 17)
        std::array<Tensor<double, 2>, 16> G_ijmn_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                G_ijmn_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("G_ijmn", nlmo_ijkl, nlmo_ijkl);
                // Jiang Eq. 17a (mi|nj)
                einsum(0.0, Indices{index::m, index::n}, &G_ijmn_list[i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::m},
                        q_io_t1_list[i_idx], Indices{index::Q, index::n}, q_io_t1_list[j_idx]);

                // Jiang Eq. 17b (me|nf) T_{ij}^{ef}
                einsum(1.0, Indices{index::m, index::n}, &G_ijmn_list[i_idx * FOUR + j_idx], 1.0, 
                        Indices{index::m, index::n, index::e, index::f}, g_menf_t, Indices{index::e, index::f}, T_ij_list[i_idx * FOUR + j_idx]);
            } // end j_idx
        } // end i_idx

        // H_{ef}^{ab} (dimensions: (e, f, a, b)) (Jiang and Matthews Eq. 19)
        Tensor<double, 4> H_efab("H_efab", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl); {
            // Jiang Eq. 19a (ae|bf)
            Tensor<double, 4> H_efab_temp("H_efab_temp", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::a, index::f, index::b}, &H_efab_temp, 1.0, Indices{index::Q, index::e, index::a}, q_vv_t1, 
                    Indices{index::Q, index::f, index::b}, q_vv_t1);
            permute(Indices{index::e, index::f, index::a, index::b}, &H_efab, Indices{index::e, index::a, index::f, index::b}, H_efab_temp);

            // Jiang Eq. 19b (me|nf)T_{mn}^{ab}
            einsum(1.0, Indices{index::e, index::f, index::a, index::b}, &H_efab, 1.0, 
                    Indices{index::m, index::n, index::e, index::f}, g_menf_t, Indices{index::m, index::n, index::a, index::b}, T_mn);
        }

        // I_{eij}^{mab} (dimensions: 16 * (m, e, a, b)) (Jiang and Matthews Eq. 21)
        std::array<Tensor<double, 4>, 16> I_eijmab_list; {
            std::array<Tensor<double, 4>, 16> I_eijmab_temp; 
            for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
                for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                    I_eijmab_temp[i_idx * FOUR + j_idx] = Tensor<double, 4>("I_eijmab", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);

                    // Jiang Eq. 21a [2(me|af) - (mf|ae)]T_{ij}^{fb}
                    einsum(0.0, Indices{index::m, index::e, index::a, index::b}, &I_eijmab_temp[i_idx * FOUR + j_idx], 1.0, 
                            Indices{index::m, index::e, index::a, index::f}, L_mfae_t, Indices{index::f, index::b}, T_ij_list[i_idx * FOUR + j_idx]);

                    // Jiang Eq. 21b -[2(me|ni) - (mi|ne)]T_{nj}^{ab}
                    einsum(1.0, Indices{index::m, index::e, index::a, index::b}, &I_eijmab_temp[i_idx * FOUR + j_idx], -1.0, 
                            Indices{index::m, index::e, index::n}, L_menj_t_list[i_idx], Indices{index::n, index::a, index::b}, T_mi_list[j_idx]);
                } // end j_idx
            } // end i_idx

            // Symmetrization P_{ij}^{ab}
            for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
                for (int j_idx = i_idx; j_idx < FOUR; ++j_idx) {
                    I_eijmab_list[i_idx * FOUR + j_idx] = I_eijmab_temp[i_idx * FOUR + j_idx];
                    Tensor<double, 4> I_ejimba_buffer = Tensor<double, 4>("I_ejimba", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::m, index::e, index::b, index::a}, &I_ejimba_buffer, 
                            Indices{index::m, index::e, index::a, index::b}, I_eijmab_temp[j_idx * FOUR + i_idx]);
                    I_eijmab_list[i_idx * FOUR + j_idx] += I_ejimba_buffer;

                    // Jiang Eq. 21c 0.25 [2(me|nf) - (mf|ne)] Z_{nij}^{fab} -> O(N^{10}) worst case (x2 for permutational symmetry)
                    einsum(1.0, Indices{index::m, index::e, index::a, index::b}, &I_eijmab_list[i_idx * FOUR + j_idx], 0.5,
                            Indices{index::m, index::e, index::n, index::f}, L_menf, Indices{index::n, index::f, index::a, index::b}, Z_mij_list[i_idx * FOUR + j_idx]);
                } // end j_idx
            } // end i_idx
        }

        // J_{iej}^{mab} (dimensions: 16 * (m, e, a, b)) (Jiang and Matthews Eq. 23)
        std::array<Tensor<double, 4>, 16> J_iejmab_list; {
            Tensor<double, 4> g_mfae_v("g_mfae_v", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::e, index::a, index::f}, &g_mfae_v, Indices{index::m, index::f, index::e, index::a}, g_mfae);
            Tensor<double, 4> g_menf_v("g_menf_v", nlmo_ijkl, nqno_ijkl, nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::e, index::n, index::f}, &g_menf_v, Indices{index::m, index::f, index::n, index::e}, g_menf);

            for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
                for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                    J_iejmab_list[i_idx * FOUR + j_idx] = Tensor<double, 4>("J_iejmab", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);

                    // Jiang Eq. 23a [(mf|ae)]T_{ij}^{fb}
                    einsum(0.0, Indices{index::m, index::e, index::a, index::b}, &J_iejmab_list[i_idx * FOUR + j_idx], 1.0, 
                            Indices{index::m, index::e, index::a, index::f}, g_mfae_v, Indices{index::f, index::b}, T_ij_list[i_idx * FOUR + j_idx]);

                    // Jiang Eq. 23b -[(mi|ne)]T_{nj}^{ab}
                    Tensor<double, 3> g_menj_u("g_menj_u", nlmo_ijkl, nqno_ijkl, nlmo_ijkl);
                    permute(Indices{index::m, index::e, index::n}, &g_menj_u, Indices{index::e, index::n, index::m}, g_menj_list[i_idx]);
                    einsum(1.0, Indices{index::m, index::e, index::a, index::b}, &J_iejmab_list[i_idx * FOUR + j_idx], -1.0,
                            Indices{index::m, index::e, index::n}, g_menj_u, Indices{index::n, index::a, index::b}, T_mi_list[j_idx]);

                    Tensor<double, 4> T_mij_t("T_mij_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::n, index::f, index::a, index::b}, &T_mij_t, Indices{index::n, index::a, index::f, index::b}, T_mij_list[i_idx * FOUR + j_idx]);

                    // Jiang Eq. 23c -0.5 (mf|ne) T_{nij}^{afb} -> O(N^{10}) worst case
                    einsum(1.0, Indices{index::m, index::e, index::a, index::b}, &J_iejmab_list[i_idx * FOUR + j_idx], -0.5,
                            Indices{index::m, index::e, index::n, index::f}, g_menf_v, Indices{index::n, index::f, index::a, index::b}, T_mij_t);
                } // end j_idx
            } // end i_idx
        } // end scope

        // K_{ijk}^{amn} (dimensions: 24 * (a, m, n)) (Jiang and Matthews Eq. 25)
        std::unordered_map<int, Tensor<double, 3>> K_ijkamn_map; {
            std::unordered_map<int, Tensor<double, 3>> K_ijkamn_temp;

            einsums::for_sequence<24UL>([&](auto perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
                int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
                int ijk_idx = i * std::pow(naocc, 2) + j * naocc + k;

                K_ijkamn_temp[perm_idx] = Tensor<double, 3>("K_ijkamn", nqno_ijkl, nlmo_ijkl, nlmo_ijkl);

                // Jiang Eq. 25a (me|nk)T_{ij}^{ae}
                einsum(0.0, Indices{index::a, index::m, index::n}, &K_ijkamn_temp[perm_idx], 1.0,
                        Indices{index::a, index::e}, T_ij_list[i_idx * FOUR + j_idx], Indices{index::e, index::m, index::n}, g_menj_list[k_idx]);

                // Jiang Eq. 25b 0.5 T_{ijk}^{aef} (me|nf)
                if (i_j_k_to_ijk_.count(ijk_idx)) {
                    int ijk = i_j_k_to_ijk_[ijk_idx];
                    auto S_ijkl_ijk = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[ijk]);
                    S_ijkl_ijk = linalg::triplet(X_qno_[ijkl], S_ijkl_ijk, X_tno_[ijk], true, false, false);

                    Tensor<double, 3> T_ijk = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[ijk], i, j, k), S_ijkl_ijk, n_tno_[ijk], n_qno_[ijkl]);
                
                    einsum(1.0, Indices{index::a, index::m, index::n}, &K_ijkamn_temp[perm_idx], 0.5, 
                            Indices{index::a, index::e, index::f}, T_ijk, Indices{index::m, index::n, index::e, index::f}, g_menf_t);
                } // end if
            });

            // (1 + P_{jk}^{mn})
            einsums::for_sequence<24UL>([&](auto perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
                int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
                int ijk_idx = i * std::pow(naocc, 2) + j * naocc + k;

                K_ijkamn_map[perm_idx] = K_ijkamn_temp[perm_idx];

                size_t twin_perm_idx = ijk_to_ijkl_perm_idx[i * std::pow(naocc, 2) + k * naocc + j];
                Tensor<double, 3> K_ikjanm_temp("K_ikjanm_temp", nqno_ijkl, nlmo_ijkl, nlmo_ijkl);
                permute(Indices{index::a, index::n, index::m}, &K_ikjanm_temp, Indices{index::a, index::m, index::n}, K_ijkamn_temp[twin_perm_idx]);

                K_ijkamn_map[perm_idx] += K_ikjanm_temp;
            });
        }

        // L_{ijk}^{abm} (dimensions: 24 * (m, a, b)) (Jiang and Matthews Eq. 27)
        std::unordered_map<int, Tensor<double, 3>> L_ijkabm_map; {
            std::unordered_map<int, Tensor<double, 3>> L_ijkabm_temp;

            einsums::for_sequence<24UL>([&](auto perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
                int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
                int ij = i_j_to_ij_[i][j];
                int ijk_idx = i * std::pow(naocc, 2) + j * naocc + k;

                Tensor<double, 3> F_mae("F_mae", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::m, index::a, index::e}, &F_mae, Indices{index::m, index::e, index::a}, F_iema_list[k_idx]);
                
                Tensor<double, 3> EF_mea_sum = E_eima_list[i_idx]; // (m, e, a) dimensions
                EF_mea_sum += F_iema_list[i_idx];
                Tensor<double, 3> EF_mae("EF_mae", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::m, index::a, index::e}, &EF_mae, Indices{index::m, index::e, index::a}, EF_mea_sum);
                    
                L_ijkabm_temp[perm_idx] = Tensor<double, 3>("L_ijkabm_temp", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                L_ijkabm_temp[perm_idx].zero();

                if (i_j_k_to_ijk_.count(ijk_idx)) {
                    int ijk = i_j_k_to_ijk_[ijk_idx];
                    auto S_ijkl_ijk = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[ijk]);
                    S_ijkl_ijk = linalg::triplet(X_qno_[ijkl], S_ijkl_ijk, X_tno_[ijk], true, false, false);

                    Tensor<double, 3> T_ijk = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[ijk], i, j, k), S_ijkl_ijk, n_tno_[ijk], n_qno_[ijkl]);

                    // Jiang Eq. 27a (mf|ae) T_{ijk}^{ebf}
                    Tensor<double, 3> T_ijk_t("T_ijk_t", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::f, index::e, index::b}, &T_ijk_t, Indices{index::e, index::b, index::f}, T_ijk);
                    einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm_temp[perm_idx], 1.0, 
                            Indices{index::m, index::a, index::f, index::e}, g_mfae_u, Indices{index::f, index::e, index::b}, T_ijk_t);
                }

                // Jiang Eq. 27b 0.5 (E_{ei}^{ma} + F_{ie}^{ma}) T_{jk}^{be}
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm_temp[perm_idx], 0.5,
                        Indices{index::m, index::a, index::e}, EF_mae, Indices{index::b, index::e}, T_ij_list[j_idx * FOUR + k_idx]);

                // Jiang Eq. 27c (F_{ke}^{ma} T_{ij}^{eb})
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm_temp[perm_idx], 1.0,
                        Indices{index::m, index::a, index::e}, F_mae, Indices{index::e, index::b}, T_ij_list[i_idx * FOUR + j_idx]);

                // Jiang Eq. 27d -0.5 G_{ki}^{mn} T_{nj}^{ab}
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm_temp[perm_idx], -0.5,
                        Indices{index::m, index::n}, G_ijmn_list[k_idx * FOUR + i_idx], Indices{index::n, index::a, index::b}, T_mi_list[j_idx]);

                // Jiang Eq. 27e 0.5 (me|nf) alpha_{nijk}^{fabe} -> O(N^{10})
                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);
                Tensor<double, 2> S_ijkl_ij_ein("S_ijkl_ij_ein", nqno_ijkl, n_pno_[ij]);
                ::memcpy(S_ijkl_ij_ein.data(), S_ijkl_ij->get_pointer(), nqno_ijkl * n_pno_[ij] * sizeof(double));

                Tensor<double, 3> L_ijkabm_cont("L_ijkabm_cont", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                L_ijkabm_cont.zero();
                int k_ij = lmopair_to_lmos_dense_[ij][k];

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    int m_ij = lmopair_to_lmos_dense_[ij][m];

                    if (m_ij == -1) continue;

                    Tensor<double, 2> L_alpha_slice = (L_alpha_list[ij])(k_ij, m_ij, All, All);
                    Tensor<double, 2> L_alpha_qno_pno("L_alpha_qno_pno", nqno_ijkl, n_pno_[ij]);
                    einsum(0.0, Indices{index::c, index::b}, &L_alpha_qno_pno, 1.0, Indices{index::c, index::a}, S_ijkl_ij_ein,
                            Indices{index::a, index::b}, L_alpha_slice);
                    Tensor<double, 2> L_alpha_qno_qno("L_alpha_qno_qno", nqno_ijkl, nqno_ijkl);
                    einsum(0.0, Indices{index::a, index::d}, &L_alpha_qno_qno, 1.0, Indices{index::a, index::b}, L_alpha_qno_pno,
                            Indices{index::d, index::b}, S_ijkl_ij_ein);

                    (L_ijkabm_cont)(m_ijkl, All, All) = L_alpha_qno_qno;
                } // end m_ijkl

                L_ijkabm_temp[perm_idx] += L_ijkabm_cont;
            });

            // (1 + P_{ij}^{ab})
            einsums::for_sequence<24UL>([&](auto perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
                int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
                int ijk_idx = i * std::pow(naocc, 2) + j * naocc + k;

                L_ijkabm_map[perm_idx] = L_ijkabm_temp[perm_idx];

                size_t twin_perm_idx = ijk_to_ijkl_perm_idx[j * std::pow(naocc, 2) + i * naocc + k];
                Tensor<double, 3> L_jikbam("L_jikbam", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::m, index::a, index::b}, &L_jikbam, Indices{index::m, index::b, index::a}, L_ijkabm_temp[twin_perm_idx]);
                L_ijkabm_map[perm_idx] += L_jikbam;
            });
        } // end k_ij

        // M_{ejk}^{abc} (dimension: 16 * (e, a, b, c)) (Jiang and Matthews Eq. 29)
        std::array<Tensor<double, 4>, 16> M_ejkabc_list; {
            std::array<Tensor<double, 4>, 16> M_ejkabc_temp;
            Tensor<double, 4> H_efab_transpose("H_efab_transpose", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::e, index::a, index::b, index::f}, &H_efab_transpose, Indices{index::e, index::f, index::a, index::b}, H_efab);

            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                int j = ijkl_list[j_idx];
                for (int k_idx = 0; k_idx < FOUR; ++k_idx) {
                    int k = ijkl_list[k_idx];
                    int jk = i_j_to_ij_[j][k];

                    // Jiang Eq. 29a 0.5 * H_{ef}^{ab} T_{jk}^{fc}
                    M_ejkabc_temp[j_idx * FOUR + k_idx] = Tensor<double, 4>("M_ejkabc_temp", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    einsum(0.0, Indices{index::e, index::a, index::b, index::c}, &M_ejkabc_temp[j_idx * FOUR + k_idx], 0.5,
                            Indices{index::e, index::a, index::b, index::f}, H_efab_transpose, Indices{index::f, index::c}, T_ij_list[j_idx * FOUR + k_idx]);

                    // Jiang Eq. 29b -0.5 (me|nf) alpha_{nmjk}^{fabc}
                    auto S_ijkl_jk = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[jk]);
                    S_ijkl_jk = linalg::triplet(X_qno_[ijkl], S_ijkl_jk, X_pno_[jk], true, false, false);
                    M_ejkabc_temp[j_idx * FOUR + k_idx] += matmul_4d(M_alpha_list[jk], S_ijkl_jk, n_pno_[jk], n_qno_[ijkl]);
                } // end k_idx
            } // j_idx

            // (1 + P_{jk}^{bc})
            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                for (int k_idx = 0; k_idx < FOUR; ++k_idx) {
                    M_ejkabc_list[j_idx * FOUR + k_idx] = M_ejkabc_temp[j_idx * FOUR + k_idx];

                    Tensor<double, 4> M_ekjacb("M_ekjacb", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::e, index::a, index::b, index::c}, &M_ekjacb, Indices{index::e, index::a, index::c, index::b}, M_ejkabc_temp[k_idx * FOUR + j_idx]);
                    M_ejkabc_list[j_idx * FOUR + k_idx] += M_ekjacb;
                }
            }
        }

        // => Form all possible R_ijkl's over unique (i, j, k, l)
        Tensor<double, 4> R_ijkl("R_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
        R_ijkl.zero();

        std::unordered_map<int, Tensor<double, 4>> R_ijkl_list;

        // Terms with (i, jkl)-type symmetry
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            auto &[j, k, l] = exclusion_list[i_idx];

            Tensor<double, 4> R_ijkl_buffer("R_ijkl_buffer", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            R_ijkl_buffer.zero();

            // Jiang and Matthews Eq. 12 (1/6 * C_{ae}) T_{ijkl}^{ebcd}
            // (permutationally adapted coefficient is +1)
            Tensor<double, 4> T_ijkl = quadruples_permuter(T_iajbkcld_[ijkl], i, j, k, l);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer, 1.0,
                    Indices{index::a, index::e}, C_ae, Indices{index::e, index::b, index::c, index::d}, T_ijkl);

            // Jiang and Matthews Eq. 12 (-1/6 D_{mi}) T_{mjkl}^{abcd} (permutationally adapted coefficient is -1)
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer, -1.0,
                    Indices{index::m, index::a, index::b, index::c, index::d}, T_nijk_exclusion_list[i_idx], Indices{index::m}, D_mi_list[i_idx]);

            // Jiang and Matthews Eq. 14 (1/12 E_{ei}^{ma} alpha_{mjkl}^{ebcd}) -> O(N^{10}) worst case (permutationally adapted coefficient is 1/2)
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer, 0.5,
                    Indices{index::m, index::e, index::a}, E_eima_list[i_idx], Indices{index::m, index::e, index::b, index::c, index::d}, alpha_nijk_exclusion_list[i_idx]);

            if (i_idx == 0) {
                R_ijkl += R_ijkl_buffer;
            } else if (i_idx == 1) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer, 1, 0, 2, 3);
            } else if (i_idx == 2) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer, 2, 1, 0, 3);
            } else {
                R_ijkl += quadruples_permuter(R_ijkl_buffer, 3, 1, 2, 0);
            }
        }

        // Terms with (ij, kl)-type symmetry
        std::array<std::tuple<int, int>, 6> ijkl_pair_idx = {std::make_tuple(0, 1), std::make_tuple(0, 2), std::make_tuple(0, 3), 
            std::make_tuple(1, 2), std::make_tuple(1, 3), std::make_tuple(2, 3)};
        std::array<std::tuple<int, int>, 6> ijkl_pair_idx_complement = {std::make_tuple(2, 3), std::make_tuple(1, 3), std::make_tuple(1, 2), 
            std::make_tuple(0, 3), std::make_tuple(0, 2), std::make_tuple(0, 1)};

        for (int ij_idx = 0; ij_idx < ijkl_pair_idx.size(); ++ij_idx) {
            auto &[i_idx, j_idx] = ijkl_pair_idx[ij_idx];
            auto &[k_idx, l_idx] = ijkl_pair_idx_complement[ij_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            int kl = i_j_to_ij_[k][l];

            Tensor<double, 4> R_ijkl_buffer_a("R_ijkl_buffer_a", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 4> R_ijkl_buffer_b("R_ijkl_buffer_b", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 4> R_ijkl_buffer_c("R_ijkl_buffer_c", n_pno_ext_[kl], n_pno_ext_[kl], n_pno_ext_[kl], n_pno_ext_[kl]);

            R_ijkl_buffer_a.zero();
            R_ijkl_buffer_b.zero();
            R_ijkl_buffer_c.zero();

            // Jiang and Matthews Eq. 8 (0.5 * A_{ej}^{ab} T_{ikl}^{ecd})
            int ikl_dense = i * naocc * naocc + k * naocc + l;
            if (i_j_k_to_ijk_.count(ikl_dense)) {
                // (j, i, k, l) contribution
                int ikl = i_j_k_to_ijk_[ikl_dense];
                auto S_ijkl_ikl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[ikl]);
                S_ijkl_ikl = linalg::triplet(X_qno_[ijkl], S_ijkl_ikl, X_tno_[ikl], true, false, false);
                auto T_ikl = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[ikl], i, k, l), 
                                                S_ijkl_ikl, n_tno_[ikl], n_qno_[ijkl]);
                einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, 
                        Indices{index::e, index::a, index::b}, A_ejab_list[j_idx], Indices{index::e, index::c, index::d}, T_ikl);
            }

            int jkl_dense = j * naocc * naocc + k * naocc + l;
            if (i_j_k_to_ijk_.count(jkl_dense)) {
                // (i, j, k, l) contribution
                int jkl = i_j_k_to_ijk_[jkl_dense];
                auto S_ijkl_jkl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[jkl]);
                S_ijkl_jkl = linalg::triplet(X_qno_[ijkl], S_ijkl_jkl, X_tno_[jkl], true, false, false);
                auto T_jkl = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[jkl], j, k, l), 
                                                S_ijkl_jkl, n_tno_[jkl], n_qno_[ijkl]);
                Tensor<double, 3> A_tilde("A_tilde", n_qno_[ijkl], n_qno_[ijkl], n_qno_[ijkl]);
                permute(Indices{index::e, index::a, index::b}, &A_tilde, Indices{index::e, index::b, index::a}, A_ejab_list[i_idx]);
                einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, 
                        Indices{index::e, index::a, index::b}, A_tilde, Indices{index::e, index::c, index::d}, T_jkl);
            }

            // Jiang and Matthews Eq. 10 (-0.5 * B_{ij}^{am} T_{mkl}^{bcd})
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, -1.0,
                    Indices{index::m, index::a}, B_ijam_list[i_idx * FOUR + j_idx], Indices{index::m, index::b, index::c, index::d}, T_mij_list[k_idx * FOUR + l_idx]);

            einsum(0.0, Indices{index::b, index::a, index::c, index::d}, &R_ijkl_buffer_b, -1.0,
                    Indices{index::m, index::b}, B_ijam_list[j_idx * FOUR + i_idx], Indices{index::m, index::a, index::c, index::d}, T_mij_list[k_idx * FOUR + l_idx]);

            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 1, 0, 2, 3);

            // Jiang and Matthews Eq. 16 -0.5 (0.5 + P_ab) F_{ie}^{ma} T_{jmkl}^{ebcd} -> O(N^{10}) worst case
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                // Order both (m, i, j, l) and (m, j, k, l) for the permutation lookup.
                Tensor<double, 4> T_mikl = T_nijk_exclusion_list_unsorted[j_idx](m_ijkl, All, All, All, All);
                Tensor<double, 4> T_mjkl = T_nijk_exclusion_list_unsorted[i_idx](m_ijkl, All, All, All, All);

                Tensor<double, 4> T_imkl = quadruples_permuter(T_mikl, i, m, k, l);
                Tensor<double, 4> T_jmkl = quadruples_permuter(T_mjkl, j, m, k, l);

                Tensor<double, 2> F_iema_slice = F_iema_list[i_idx](m_ijkl, All, All);
                Tensor<double, 2> F_jemb_slice = F_iema_list[j_idx](m_ijkl, All, All);

                einsum(0.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_b, -0.5, Indices{index::e, index::a}, 
                        F_iema_slice, Indices{index::e, index::b, index::c, index::d}, T_jmkl);
                R_ijkl_buffer_a += R_ijkl_buffer_b;
                R_ijkl_buffer_b = quadruples_permuter(R_ijkl_buffer_b, 1, 0, 2, 3);
                R_ijkl_buffer_b *= 2.0;
                R_ijkl_buffer_a += R_ijkl_buffer_b;

                einsum(0.0, Indices{index::b, index::a, index::c, index::d}, &R_ijkl_buffer_b, -0.5, Indices{index::e, index::b}, 
                        F_jemb_slice, Indices{index::e, index::a, index::c, index::d}, T_imkl);
                R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 1, 0, 2, 3);
                R_ijkl_buffer_b *= 2.0;
                R_ijkl_buffer_a += R_ijkl_buffer_b;
            }
            
            // Jiang and Matthews Eq. 18 0.25 * G_{ij}^{mn} T_{mnkl}^{abcd} -> O(N^{10}) worst case
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                
                for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                    int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];
                    int mn = i_j_to_ij_[m][n];
                    int mnkl_idx = m * std::pow(naocc, 3) + n * std::pow(naocc, 2) + k * naocc + l;
                    int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;
                    if (mn == -1 || mnkl == -1) continue;
                    
                    Tensor<double, 4> T_mnkl("T_mnkl", n_pno_ext_[kl], n_pno_ext_[kl], n_pno_ext_[kl], n_pno_ext_[kl]);
                    int lk = ij_to_ji_[kl], nm = ij_to_ji_[mn];

                    if (m > n && k > l) {
                        permute(Indices{index::b, index::a, index::d, index::c}, &T_mnkl, Indices{index::a, index::b, index::c, index::d}, T_mnkl_list_[lk][nm]);
                    } else if (m > n) {
                        permute(Indices{index::b, index::a, index::c, index::d}, &T_mnkl, Indices{index::a, index::b, index::c, index::d}, T_mnkl_list_[kl][nm]);
                    } else if (k > l) {
                        permute(Indices{index::a, index::b, index::d, index::c}, &T_mnkl, Indices{index::a, index::b, index::c, index::d}, T_mnkl_list_[lk][mn]);
                    } else {
                        T_mnkl = T_mnkl_list_[kl][mn];
                    }

                    size_t length = std::pow(n_pno_ext_[kl], 4);
                    C_DAXPY(length, (G_ijmn_list[i_idx * FOUR + j_idx])(m_ijkl, n_ijkl), T_mnkl.data(), 1, R_ijkl_buffer_c.data(), 1);
                } // end n_ijkl
            } // end m_ijkl
            // Flush G contributions
            auto S_ijkl_kl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_ext_[kl]);
            S_ijkl_kl = linalg::triplet(X_qno_[ijkl], S_ijkl_kl, X_pno_ext_[kl], true, false, false);
            R_ijkl_buffer_a += matmul_4d(R_ijkl_buffer_c, S_ijkl_kl, n_pno_ext_[kl], n_qno_[ijkl]);

            // Jiang and Matthews Eq. 20 0.25 * H_{ef}^{ab} T_{ijkl}^{efcd} -> O(N^{10}) worst case
            Tensor<double, 4> T_ijkl = quadruples_permuter(T_iajbkcld_[ijkl], i, j, k, l);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, Indices{index::e, index::f, index::a, index::b}, H_efab,
                    Indices{index::e, index::f, index::c, index::d}, T_ijkl);

            // Jiang and Matthews Eq. 22 0.125 * I_{eij}^{mab} Z_{mkl}^{ecd} -> O(N^{10}) worst case
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 0.5, Indices{index::m, index::e, index::a, index::b}, 
                    I_eijmab_list[i_idx * FOUR + j_idx], Indices{index::m, index::e, index::c, index::d}, Z_mij_list[k_idx * FOUR + l_idx]);

            // Jiang and Matthews Eq. 26 0.5 K_{ijk}^{amn} T_{mnl}^{bcd} -> O(N^{10}) worst case
            // (becomes) T_{mnk}^{abc} K_{lij}^{dmn}
            size_t lij_perm_idx = ijk_to_ijkl_perm_idx[l * std::pow(naocc, 2) + i * naocc + j];
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, Indices{index::m, index::n, index::a, index::b, index::c}, 
                    T_mni_list[k_idx], Indices{index::d, index::m, index::n}, K_ijkamn_map[lij_perm_idx]);

            size_t kij_perm_idx = ijk_to_ijkl_perm_idx[k * std::pow(naocc, 2) + i * naocc + j];
            einsum(0.0, Indices{index::a, index::b, index::d, index::c}, &R_ijkl_buffer_b, 1.0, Indices{index::m, index::n, index::a, index::b, index::d},
                    T_mni_list[l_idx], Indices{index::c, index::m, index::n}, K_ijkamn_map[kij_perm_idx]);
            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 0, 1, 3, 2);

            // Jiang and Matthews Eq. 28 -0.5 L_{ijk}^{abm} T_{ml}^{cd}
            size_t ijk_perm_idx = ijk_to_ijkl_perm_idx[i * std::pow(naocc, 2) + j * naocc + k];
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, -1.0, Indices{index::m, index::a, index::b}, 
                    L_ijkabm_map[ijk_perm_idx], Indices{index::m, index::c, index::d}, T_mi_list[l_idx]);
            
            size_t ijl_perm_idx = ijk_to_ijkl_perm_idx[i * std::pow(naocc, 2) + j * naocc + l];
            einsum(0.0, Indices{index::a, index::b, index::d, index::c}, &R_ijkl_buffer_b, -1.0, Indices{index::m, index::a, index::b}, 
                    L_ijkabm_map[ijl_perm_idx], Indices{index::m, index::d, index::c}, T_mi_list[k_idx]);
            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 0, 1, 3, 2);

            // Jiang and Matthews Eq. 30 0.5 M_{ejk}^{abc} T_{il}^{ed}
            // (becomes) 0.5 T_{ji}^{ea} M_{ekl}^{bcd}
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, Indices{index::e, index::a}, 
                    T_ij_list[j_idx * FOUR + i_idx], Indices{index::e, index::b, index::c, index::d}, M_ejkabc_list[k_idx * FOUR + l_idx]);
            einsum(0.0, Indices{index::b, index::a, index::c, index::d}, &R_ijkl_buffer_b, 1.0, Indices{index::e, index::b}, 
                    T_ij_list[i_idx * FOUR + j_idx], Indices{index::e, index::a, index::c, index::d}, M_ejkabc_list[k_idx * FOUR + l_idx]);
            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 1, 0, 2, 3);

            if (ij_idx == 0) {
                R_ijkl += R_ijkl_buffer_a; // (a, b, c, d)
            } else if (ij_idx == 1) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer_a, 0, 2, 1, 3); // (a, c, b, d)
            } else if (ij_idx == 2) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer_a, 0, 2, 3, 1); // (a, d, b, c)
            } else if (ij_idx == 3) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer_a, 2, 0, 1, 3); // (b, c, a, d)
            } else if (ij_idx == 4) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer_a, 2, 0, 3, 1); // (b, d, a, c)
            } else if (ij_idx == 5) {
                R_ijkl += quadruples_permuter(R_ijkl_buffer_a, 2, 3, 0, 1); // (c, d, a, b)
            }
        }

        // Loop over all possible quadruplet permutations
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            int ijkl_idx = i * std::pow(naocc, 3) + j * std::pow(naocc, 2) + k * naocc + l;

            // Accumulate the permuted Eq. 24 contributions in the quadruples residual.
            if (!R_ijkl_list.count(ijkl_idx)) {
                Tensor<double, 4> R_ijkl_perm("R_ijkl_perm", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                R_ijkl_perm.zero();

                // => Buffer used to keep track of contributions <= //
                Tensor<double, 4> R_ijkl_buffer_a("R_ijkl_buffer_a", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                Tensor<double, 4> R_ijkl_buffer_b("R_ijkl_buffer_b", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);

                // Jiang and Matthews Eq. 24 -(0.5 + P_{ac}) J_{iej}^{mab} T_{mkl}^{ced} -> O(N^{10}) worst case
                Tensor<double, 4> T_mkl_t = T_mij_list[k_idx * FOUR + l_idx];
                permute(Indices{index::m, index::e, index::c, index::d}, &T_mkl_t, Indices{index::m, index::c, index::e, index::d}, T_mij_list[k_idx * FOUR + l_idx]);

                einsum(0.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, Indices{index::m, index::e, index::a, index::b},
                        J_iejmab_list[i_idx * FOUR + j_idx], Indices{index::m, index::e, index::c, index::d}, T_mkl_t);
                permute(Indices{index::c, index::b, index::a, index::d}, &R_ijkl_buffer_b, Indices{index::a, index::b, index::c, index::d}, R_ijkl_buffer_a);

                R_ijkl_buffer_a *= -0.5;
                R_ijkl_perm += R_ijkl_buffer_a;

                R_ijkl_buffer_b *= -1.0;
                R_ijkl_perm += R_ijkl_buffer_b;

                // Add permutation to buffer
                R_ijkl_list[ijkl_idx] = R_ijkl_perm;
                permute(Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, std::get<perm_idx>(einsum_indices), R_ijkl_list[ijkl_idx]);
                R_ijkl += R_ijkl_buffer_a;
            } else {
                Tensor<double, 4> R_ijkl_buffer_a("R_ijkl_buffer_a", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, std::get<perm_idx>(einsum_indices), R_ijkl_list[ijkl_idx]);
                R_ijkl += R_ijkl_buffer_a;
            }
        }); // end if

        if (disk_ints_quads_) {
            q_ov_ijkl_[ijkl] = Tensor<double, 3>(q_ov_ijkl_[ijkl].name(), 0, 0, 0);
            q_vv_ijkl_[ijkl] = Tensor<double, 3>(q_vv_ijkl_[ijkl].name(), 0, 0, 0);
        }

        R_iajbkcld[ijkl] = R_ijkl;

    } // end ijkl
}

void DLPNOCCSDTQ::lccsdtq_iterations() {

    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    int n_lmo_triplets = ijk_to_i_j_k_.size();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    const bool extrapolate_t4 = options_.get_bool("EXTRAPOLATE_T4");

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
    std::vector<Tensor<double, 4>> R_iajbkcld(n_lmo_quadruplets);

    // Psi4 Matrix copies are required by the current DIIS interface.
    std::vector<SharedMatrix> R_iajbkcld_psi;
    std::vector<SharedMatrix> T_iajbkcld_psi;
    if (extrapolate_t4) {
        R_iajbkcld_psi.resize(n_lmo_quadruplets);
        T_iajbkcld_psi.resize(n_lmo_quadruplets);
    }

    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        R_ia[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
    }

    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        R_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        Rn_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
    }

    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        R_iajbkc[ijk] = std::make_shared<Matrix>(n_tno_[ijk], n_tno_[ijk] * n_tno_[ijk]);
    }

    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
        R_iajbkcld[ijkl] = Tensor<double, 4>("R_iajbkcld", n_qno_[ijkl], n_qno_[ijkl], n_qno_[ijkl], n_qno_[ijkl]);
        if (extrapolate_t4) {
            R_iajbkcld_psi[ijkl] = std::make_shared<Matrix>(n_qno_[ijkl] * n_qno_[ijkl], n_qno_[ijkl] * n_qno_[ijkl]);
            T_iajbkcld_psi[ijkl] = std::make_shared<Matrix>(n_qno_[ijkl] * n_qno_[ijkl], n_qno_[ijkl] * n_qno_[ijkl]);
        }
    }

    std::vector<std::vector<SharedMatrix>> R_ia_buffer(nthreads);
    std::vector<std::vector<SharedMatrix>> R_iajb_buffer(nthreads);

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
    T_n_ijkl_.resize(n_lmo_quadruplets);

    // LCCSDTQ iterations

    outfile->Printf("\n  ==> Local CCSDTQ <==\n\n");
    
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                        Corr. Energy    Delta E    RMS R1     RMS R2     RMS R3     RMS R4     Time (s)\n");

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, e_weak = 0.0, r_curr1 = 1.0, r_curr2 = 1.0, r_curr3 = 1.0, r_curr4 = 1.0;
    bool e_converged = false, r_converged = false;
    const int n_microiterations = options_.get_int("DLPNO_QUADS_MICROITERATIONS");
    
    DIISManager diis = DIISManager(options_.get_int("DIIS_MAX_VECS"), "LCCSDTQ DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::OnDisk);

    while (!(e_converged && r_converged)) {
        // RMS of residual per LMO orbital, for assessing convergence
        std::vector<double> R_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);
        // RMS of residual per LMO triplet, for assessing convergence
        std::vector<double> R_iajbkc_rms(n_lmo_triplets, 0.0);
        // RMS of residual per LMO quadruplet, for assessing convergence
        std::vector<double> R_iajbkcld_rms(n_lmo_quadruplets, 0.0);

        std::time_t time_start = std::time(nullptr);

        for (int miter = 0; miter < n_microiterations; ++miter) {

            // Create T_n_ij, T_n_ijk, and T_n_ijkl intermediates
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

    #pragma omp parallel for schedule(dynamic, 1)
            for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                int ijkl = sorted_quadruplets_[ijkl_sorted];
                auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
                int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();

                T_n_ijkl_[ijkl] = Tensor<double, 2>("T_n_ijkl", nlmo_ijkl, n_qno_[ijkl]);
                
                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    int mm = i_j_to_ij_[m][m];

                    auto S_ijkl_mm = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[mm]);
                    S_ijkl_mm = linalg::triplet(X_qno_[ijkl], S_ijkl_mm, X_pno_[mm], true, false, false);

                    auto T_m_temp = linalg::doublet(S_ijkl_mm, T_ia_[m]);

                    for (int a_ijkl = 0; a_ijkl < n_qno_[ijkl]; ++a_ijkl) {
                        (T_n_ijkl_[ijkl])(m_ijkl, a_ijkl) = (*T_m_temp)(a_ijkl, 0);
                    }
                } // end m_ijkl
            } // end ijkl
            
            // Create T_iajbkc_clone intermediate
    #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];
                auto &[i, j, k] = ijk_to_i_j_k_[ijk];

                T_iajbkc_clone_[ijk] = Tensor<double, 3>("T_ijk", n_tno_[ijk], n_tno_[ijk], n_tno_[ijk]);
                ::memcpy(T_iajbkc_clone_[ijk].data(), T_iajbkc_[ijk]->get_pointer(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));

                U_iajbkc_[ijk] = triples_spin_summation(T_iajbkc_clone_[ijk]);
            } // end ijk

            // T1-dress integrals and Fock matrices
            t1_ints();
            t1_fock();

            // spin adapt and then de-adapt triples amplitudes
        #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];
                auto &[i, j, k] = ijk_to_i_j_k_[ijk];
                
                T_iajbkc_clone_[ijk] = triples_spin_desummation(triples_spin_summation(T_iajbkc_clone_[ijk]));
                U_iajbkc_[ijk] = triples_spin_summation(T_iajbkc_clone_[ijk]);
                ::memcpy(T_iajbkc_[ijk]->get_pointer(), T_iajbkc_clone_[ijk].data(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            }

            // spin-adapt and then de-adapt quadruples amplitudes
        #pragma omp parallel for schedule(dynamic, 1)
            for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                int ijkl = sorted_quadruplets_[ijkl_sorted];
                auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

                T_iajbkcld_[ijkl] = quadruples_spin_desummation(quadruples_spin_summation(T_iajbkcld_[ijkl]));
            }

            // compute quadruples amplitude
            timer_on("DLPNO-CCSDTQ : R_iajbkcld");
            if (miter == 0) {
                form_T_mnkl();
                compute_R_iajbkcld(R_iajbkcld);

            #pragma omp parallel for schedule(dynamic, 1)
                for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                    int ijkl = sorted_quadruplets_[ijkl_sorted];
                    auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

                    R_iajbkcld[ijkl] = quadruples_spin_desummation(quadruples_spin_summation(R_iajbkcld[ijkl]));
                }

                // Update quadruples amplitude
                r_curr4 = 0.0;
        #pragma omp parallel for schedule(dynamic, 1) reduction(+ : r_curr4)
                for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                    int ijkl = sorted_quadruplets_[ijkl_sorted];
                    auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

                    Tensor<double, 4> zero_tensor("zero", n_qno_[ijkl], n_qno_[ijkl], n_qno_[ijkl], n_qno_[ijkl]);
                    zero_tensor.zero();

                    double alpha = (fabs(rmsd(R_iajbkcld[ijkl], zero_tensor)) > fabs(R_iajbkcld_rms[ijkl])) ? damping_ratio_quads_ : 0.0;

                    for (int a = 0; a < n_qno_[ijkl]; ++a) {
                        for (int b = 0; b < n_qno_[ijkl]; ++b) {
                            for (int c = 0; c < n_qno_[ijkl]; ++c) {
                                for (int d = 0; d < n_qno_[ijkl]; ++d) {
                                    (T_iajbkcld_[ijkl])(a, b, c, d) -= (1.0 - alpha) * (R_iajbkcld[ijkl])(a, b, c, d) / 
                                        ((*e_qno_[ijkl])(a) + (*e_qno_[ijkl])(b) + (*e_qno_[ijkl])(c) + (*e_qno_[ijkl])(d) - (*F_lmo_)(i,i) 
                                        - (*F_lmo_)(j,j) - (*F_lmo_)(k,k) - (*F_lmo_)(l,l));
                                } // end d
                            } // end c
                        } // end b
                    } // end a

                    // Done for DIIS
                    if (extrapolate_t4) {
                        ::memcpy(R_iajbkcld_psi[ijkl]->get_pointer(), R_iajbkcld[ijkl].data(), n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * sizeof(double));
                        ::memcpy(T_iajbkcld_psi[ijkl]->get_pointer(), T_iajbkcld_[ijkl].data(), n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * sizeof(double));
                        R_iajbkcld_rms[ijkl] = R_iajbkcld_psi[ijkl]->rms();
                    } else {
                        R_iajbkcld_rms[ijkl] = rmsd(R_iajbkcld[ijkl], zero_tensor);
                    } // end else
                    r_curr4 += R_iajbkcld_rms[ijkl] * R_iajbkcld_rms[ijkl];
                }
                r_curr4 = std::sqrt(r_curr4 / n_lmo_quadruplets);
            }
            timer_off("DLPNO-CCSDTQ : R_iajbkcld");

            // compute triples amplitude
            timer_on("DLPNO-CCSDTQ : R_iajbkc");
            if (miter == 0) {
                // form_T_mnk();
                compute_R_iajbkc_quads(R_iajbkc);

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
                    double alpha = (fabs(R_iajbkc[ijk]->rms()) > fabs(R_iajbkc_rms[ijk])) ? damping_ratio_quads_ : 0.0;

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
            }
            timer_off("DLPNO-CCSDTQ : R_iajbkc");

            // compute doubles amplitude
            timer_on("DLPNO-CCSDTQ : R_iajb");
            compute_R_iajb_quads(R_iajb, Rn_iajb, R_iajb_buffer);

            // Update doubles amplitude
            r_curr2 = 0.0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : r_curr2)
            for (int ij = 0; ij < n_lmo_pairs; ++ij) {
                auto &[i, j] = ij_to_i_j_[ij];
                double alpha = (fabs(R_iajb[ij]->rms()) > fabs(R_iajb_rms[ij])) ? damping_ratio_quads_ : 0.0;

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
            timer_off("DLPNO-CCSDTQ : R_iajb");

            // compute singles amplitude
            timer_on("DLPNO-CCSDTQ : R_ia");
            compute_R_ia_triples(R_ia, R_ia_buffer);

            // Update singles amplitude
            r_curr1 = 0.0;
    #pragma omp parallel for reduction(+ : r_curr1)
            for (int i = 0; i < naocc; ++i) {
                int ii = i_j_to_ij_[i][i];
                double alpha = (fabs(R_ia[i]->rms()) > fabs(R_ia_rms[i])) ? damping_ratio_quads_ : 0.0;

                for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                    (*T_ia_[i])(a_ii, 0) -= (1.0 - alpha) * (*R_ia[i])(a_ii, 0) / (e_pno_[ii]->get(a_ii) - F_lmo_->get(i,i));
                }
                R_ia_rms[i] = R_ia[i]->rms();
                r_curr1 += R_ia_rms[i] * R_ia_rms[i];
            }
            r_curr1 = std::sqrt(r_curr1 / naocc);
            timer_off("DLPNO-CCSDTQ : R_ia");

        } // end miter

        // DIIS Extrapolation
        size_t nelements = T_ia_.size() + T_iajb_.size() + T_iajbkc_.size();
        if (extrapolate_t4) nelements += T_iajbkcld_psi.size();

        std::vector<SharedMatrix> T_vecs;
        T_vecs.reserve(nelements);
        T_vecs.insert(T_vecs.end(), T_ia_.begin(), T_ia_.end());
        T_vecs.insert(T_vecs.end(), T_iajb_.begin(), T_iajb_.end());
        T_vecs.insert(T_vecs.end(), T_iajbkc_.begin(), T_iajbkc_.end());
        if (extrapolate_t4) T_vecs.insert(T_vecs.end(), T_iajbkcld_psi.begin(), T_iajbkcld_psi.end());

        std::vector<SharedMatrix> R_vecs;
        R_vecs.reserve(nelements);
        R_vecs.insert(R_vecs.end(), R_ia.begin(), R_ia.end());
        R_vecs.insert(R_vecs.end(), R_iajb.begin(), R_iajb.end());
        R_vecs.insert(R_vecs.end(), R_iajbkc.begin(), R_iajbkc.end());
        if (extrapolate_t4) R_vecs.insert(R_vecs.end(), R_iajbkcld_psi.begin(), R_iajbkcld_psi.end());

        auto T_vecs_flat = flatten_mats(T_vecs);
        auto R_vecs_flat = flatten_mats(R_vecs);

        if (iteration == 1) {
            diis.set_error_vector_size(R_vecs_flat);
            diis.set_vector_size(T_vecs_flat);
        }

        diis.add_entry(R_vecs_flat.get(), T_vecs_flat.get());
        diis.extrapolate(T_vecs_flat.get());

        copy_flat_mats(T_vecs_flat, T_vecs);

        // Copy data from Psi4 matrices to Einsums tensors.
        if (extrapolate_t4) {
#pragma omp parallel for schedule(dynamic, 1)
            for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                int ijkl = sorted_quadruplets_[ijkl_sorted];
                ::memcpy(T_iajbkcld_[ijkl].data(), T_iajbkcld_psi[ijkl]->get_pointer(), n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * n_qno_[ijkl] * sizeof(double));
            }
        }

        // evaluate energy and convergence
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
        r_converged &= fabs(r_curr4) < options_.get_double("R_CONVERGENCE");
        e_converged = fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE");

        e_lccsdtq_ = e_curr - e_weak;
        de_weak_ = e_weak;

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @LCCSDTQ iter %3d: %16.12f %10.3e %10.3e %10.3e %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, r_curr1, r_curr2, r_curr3, r_curr4, (int)time_stop - (int)time_start);

        ++iteration;

        if (iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }
}

double DLPNOCCSDTQ::compute_energy() {
    // Run DLPNO-CCSDT(Q) as initial step
    double e_dlpno_ccsdt_q = DLPNOCCSDT_Q::compute_energy();

    timer_on("DLPNO-CCSDTQ");

    einsums::profile::initialize();

    print_header();

    if (write_qia_pno_) {
        psio_->open(PSIF_DLPNO_QIA_PNO, PSIO_OPEN_OLD);
    }

    if (write_qab_pno_) {
        psio_->open(PSIF_DLPNO_QAB_PNO, PSIO_OPEN_OLD);
    }
    
    disk_ints_quads_ = options_.get_bool("DLPNO_CCSDTQ_DISK_INTS");
    damping_ratio_quads_ = options_.get_double("QUADRUPLES_DAMPING_RATIO");

    if (disk_ints_) {
        psio_->open(PSIF_DLPNO_QIA_TNO, PSIO_OPEN_OLD);
        psio_->open(PSIF_DLPNO_QAB_TNO, PSIO_OPEN_OLD);
    }

    if (disk_ints_quads_) {
        psio_->open(PSIF_DLPNO_QIA_QNO, PSIO_OPEN_NEW);
        psio_->open(PSIF_DLPNO_QAB_QNO, PSIO_OPEN_NEW);
    }
    // Compute extended PNOs for T4 "lasagna" terms
    double xpno_tolerance = options_.get_double("T_CUT_XPNO");
    xpno_transform(xpno_tolerance);

    timer_on("DLPNO-CCSDTQ : Estimate Memory");
    estimate_memory();
    timer_off("DLPNO-CCSDTQ : Estimate Memory");

    timer_on("DLPNO-CCSDTQ : Compute Integrals");
    compute_integrals();
    timer_off("DLPNO-CCSDTQ : Compute Integrals");

    // Compute DLPNO-CCSDTQ energy
    timer_on("DLPNO-CCSDTQ : LCCSDTQ iterations");
    lccsdtq_iterations();
    timer_off("DLPNO-CCSDTQ : LCCSDTQ iterations");

    // Close the PSIO files used for triples and quadruples intermediates.
    if (write_qia_pno_) {
        psio_->close(PSIF_DLPNO_QIA_PNO, 0);
    }

    if (write_qab_pno_) {
        psio_->close(PSIF_DLPNO_QAB_PNO, 0);
    }

    if (disk_ints_) {
        psio_->close(PSIF_DLPNO_QIA_TNO, 0);
        psio_->close(PSIF_DLPNO_QAB_TNO, 0);
    }

    if (disk_ints_quads_) {
        psio_->close(PSIF_DLPNO_QIA_QNO, 0);
        psio_->close(PSIF_DLPNO_QAB_QNO, 0);
    }

    einsums::profile::report("einsums.time", false);
    einsums::profile::finalize();

    timer_off("DLPNO-CCSDTQ");

    double e_scf = variables_["SCF TOTAL ENERGY"];
    double e_ccsdtq_corr = e_lccsdtq_ + de_lccsdt_q_screened_ + de_lccsd_t_screened_ + de_weak_ +
                           de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
    double e_ccsdtq_total = e_scf + e_ccsdtq_corr;

    set_scalar_variable("CCSDTQ CORRELATION ENERGY", e_ccsdtq_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsdtq_corr);
    set_scalar_variable("CCSDTQ TOTAL ENERGY", e_ccsdtq_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsdtq_total);

    print_results();

    return e_ccsdtq_total;
}

void DLPNOCCSDTQ::print_results() {

    int naocc = i_j_to_ij_.size();
    double t1diag = 0.0;
#pragma omp parallel for reduction(+ : t1diag)
    for (int i = 0; i < naocc; ++i) {
        t1diag += T_ia_[i]->vector_dot(T_ia_[i]);
    }
    t1diag = std::sqrt(t1diag / (2.0 * naocc));
    outfile->Printf("\n  T1 Diagnostic: %8.8f \n", t1diag);
    set_scalar_variable("CC T1 DIAGNOSTIC", t1diag);

    double e_total = e_lccsdtq_ + de_lccsdt_q_screened_ + de_lccsd_t_screened_ + de_weak_ +
                     de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;

    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSDTQ Correlation Energy: %16.12f \n", e_total);
    outfile->Printf("    LCCSDTQ Correlation Energy:          %16.12f \n", e_lccsdtq_);
    outfile->Printf("    Weak Pair Contribution:              %16.12f \n", de_weak_);
    outfile->Printf("    Eliminated Pair MP2 Correction:      %16.12f \n", de_lmp2_eliminated_);
    outfile->Printf("    Dipole Pair Correction:              %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:           %16.12f \n", de_pno_total_);
    outfile->Printf("    Screened Triplets Contribution:      %16.12f \n", de_lccsd_t_screened_);
    outfile->Printf("    Screened Quadruples Correction:      %16.12f \n", de_lccsdt_q_screened_);
    outfile->Printf("\n\n  @Total DLPNO-CCSDTQ Energy: %16.12f \n", variables_["SCF TOTAL ENERGY"] + e_total);
}

} // namespace dlpno
} // namespace psi
