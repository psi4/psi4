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
#include <array>
#include <cstring>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

using einsums::All;
using einsums::Indices;
using einsums::Tensor;
using einsums::index::a;
using einsums::index::b;
using einsums::index::c;
using einsums::index::d;
using einsums::tensor_algebra::einsum;
using einsums::tensor_algebra::permute;
namespace index = einsums::index;
namespace linear_algebra = einsums::linear_algebra;

namespace {

std::string quadruples_record_name(const std::string& label, int ijkl, int component = -1) {
    std::string name = label + " " + std::to_string(ijkl);
    if (component >= 0) name += " " + std::to_string(component);
    return name;
}

size_t quadruplet_key(int i, int j, int k, int l, size_t nocc) {
    return (triplet_key(i, j, k, nocc) * nocc) + l;
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
    const bool full_quadruples_follow = algorithm_ == DLPNOMethod::CCSDTQ;
    const bool q0_only = !full_quadruples_follow && options_.get_bool("Q0_APPROXIMATION");
    const double t_cut_qno_full = options_.get_double("T_CUT_QNO_FULL");
    const double t_cut_qno_strong =
        full_quadruples_follow ? t_cut_qno_full : options_.get_double("T_CUT_QNO_STRONG");
    const double t_cut_qno_weak =
        full_quadruples_follow ? t_cut_qno_full : options_.get_double("T_CUT_QNO_WEAK");

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
    outfile->Printf("    T_CUT_QNO_STRONG                = %6.3e \n", t_cut_qno_strong);
    outfile->Printf("    T_CUT_QNO_WEAK                  = %6.3e \n", t_cut_qno_weak);
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
    outfile->Printf("    E_CONVERGENCE                   = %6.3e \n",
                    options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE                   = %6.3e \n",
                    options_.get_double("R_CONVERGENCE"));
    outfile->Printf("    DLPNO_MAXITER                   = %6d   \n",
                    options_.get_int("DLPNO_MAXITER"));
    outfile->Printf("    Q0_APPROXIMATION                = %6s   \n",
                    q0_only ? "TRUE" : "FALSE");
    outfile->Printf("    DLPNO_TOGGLE_MEMORY             = %6s   \n",
                    toggle_memory_ ? "TRUE" : "FALSE");
    outfile->Printf("    WRITE_QUADRUPLES_INTERMEDIATES  = %6s   \n\n",
                    write_quadruples_intermediates_ ? "TRUE" : "FALSE");
    if (full_quadruples_follow &&
        (options_["T_CUT_QNO_STRONG"].has_changed() || options_["T_CUT_QNO_WEAK"].has_changed())) {
        outfile->Printf(
            "    WARNING: T_CUT_QNO_STRONG/T_CUT_QNO_WEAK apply only when DLPNO-CCSDT(Q)\n"
            "             is the final method. Full quadruples were requested, so both iterative-(Q)\n"
            "             cutoffs are overridden by T_CUT_QNO_FULL = %6.3e.\n\n",
            t_cut_qno_full);
    }
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
    }
    for (size_t component = 0; component < intermediates.K_iajb.size(); ++component) {
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
        load_quads_record(psio_, quadruples_record_name("K_iajm", ijkl, component), nqno, nlmo,
                          intermediates.K_iajm[component].data());
    }
    for (size_t component = 0; component < intermediates.K_iajb.size(); ++component) {
        intermediates.K_iajb[component] = Tensor<double, 2>("K_iajb", nqno, nqno);
        intermediates.U_iajb[component] = Tensor<double, 2>("U_iajb", nqno, nqno);
        load_quads_record(psio_, quadruples_record_name("K_iajb", ijkl, component), nqno, nqno,
                          intermediates.K_iajb[component].data());
        load_quads_record(psio_, quadruples_record_name("U_iajb", ijkl, component), nqno, nqno,
                          intermediates.U_iajb[component].data());
    }
    return intermediates;
}

Tensor<double, 4> DLPNOCCSDT_Q::matmul_4d(const Tensor<double, 4>& A,
                                           const SharedMatrix& X, int dim_old,
                                           int dim_new) {
    return matmul_4d_permuted(A, X, dim_old, dim_new, 0, 1, 2, 3);
}

Tensor<double, 4> DLPNOCCSDT_Q::matmul_4d_permuted(
    const Tensor<double, 4>& A, const SharedMatrix& X, int dim_old, int dim_new,
    int i, int j, int k, int l) {
    // Transform A'[i,j,k,l] = A[I,J,K,L] X[i,I] X[j,J] X[k,K] X[l,L].
    // Before every contraction, the old index is explicitly moved to the
    // contiguous final axis. Consequently each einsum is a matrix multiply.
    // Only the current and next rank-four buffers coexist.
    Tensor<double, 2> Xview("Xview", dim_new, dim_old);
    ::memcpy(Xview.data(), X->get_pointer(),
             static_cast<size_t>(dim_new) * dim_old * sizeof(double));

    Tensor<double, 4> work("rank4_transform_1", dim_old, dim_old, dim_old, dim_new);
    einsum(0.0, Indices{index::I, index::J, index::K, index::l}, &work,
           1.0, Indices{index::I, index::J, index::K, index::L}, A,
           Indices{index::l, index::L}, Xview);

    {
        Tensor<double, 4> contiguous("rank4_transform_perm_2", dim_old, dim_old,
                                      dim_new, dim_old);
        permute(Indices{index::I, index::J, index::l, index::K}, &contiguous,
                Indices{index::I, index::J, index::K, index::l}, work);
        work = Tensor<double, 4>("released", 0, 0, 0, 0);
        Tensor<double, 4> next("rank4_transform_2", dim_old, dim_old, dim_new,
                                dim_new);
        einsum(0.0, Indices{index::I, index::J, index::l, index::k}, &next,
               1.0, Indices{index::I, index::J, index::l, index::K}, contiguous,
               Indices{index::k, index::K}, Xview);
        work = std::move(next);
    }

    {
        Tensor<double, 4> contiguous("rank4_transform_perm_3", dim_old, dim_new,
                                      dim_new, dim_old);
        permute(Indices{index::I, index::l, index::k, index::J}, &contiguous,
                Indices{index::I, index::J, index::l, index::k}, work);
        work = Tensor<double, 4>("released", 0, 0, 0, 0);
        Tensor<double, 4> next("rank4_transform_3", dim_old, dim_new, dim_new,
                                dim_new);
        einsum(0.0, Indices{index::I, index::l, index::k, index::j}, &next,
               1.0, Indices{index::I, index::l, index::k, index::J}, contiguous,
               Indices{index::j, index::J}, Xview);
        work = std::move(next);
    }

    {
        Tensor<double, 4> contiguous("rank4_transform_perm_4", dim_new, dim_new,
                                      dim_new, dim_old);
        permute(Indices{index::l, index::k, index::j, index::I}, &contiguous,
                Indices{index::I, index::l, index::k, index::j}, work);
        work = Tensor<double, 4>("released", 0, 0, 0, 0);
        Tensor<double, 4> next("rank4_transform_4", dim_new, dim_new, dim_new,
                                dim_new);
        einsum(0.0, Indices{index::l, index::k, index::j, index::i}, &next,
               1.0, Indices{index::l, index::k, index::j, index::I}, contiguous,
               Indices{index::i, index::I}, Xview);
        work = std::move(next);
    }

    Tensor<double, 4> result("rank4_transform_result", dim_new, dim_new, dim_new,
                              dim_new);
    permute(Indices{index::i, index::j, index::k, index::l}, &result,
            Indices{index::l, index::k, index::j, index::i}, work);
    work = Tensor<double, 4>("released", 0, 0, 0, 0);

    // Transforming before the occupied-column permutation makes every
    // contraction a GEMM regardless of the requested order.  Only the smaller
    // dim_new^4 result is permuted; the old-space amplitude is never copied.
    // With four distinct occupied indices, the only identity ordering is the
    // strictly increasing one. Do not use a stable-sort identity shortcut
    // here: for repeated indices, quadruples_permuter() deliberately applies
    // its historical first-match tie convention (for example, an i == j
    // column swap). Skipping that permutation changes iterative-(Q) couplings
    // for repeated-index quadruplets.
    if (i < j && j < k && k < l) return result;
    return quadruples_permuter(result, i, j, k, l);
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

    // Post-CCSD(T) methods collapse T_CUT_PAIRS_MP2 onto T_CUT_PAIRS, so every
    // pair retained after MP2 prescreening must be a strong pair. Keep this
    // invariant explicit rather than maintaining an unreachable weak-pair path.
    if (!ij_to_i_j_weak_.empty()) {
        throw PSIEXCEPTION(
            "Post-CCSD(T) quadruplet construction requires all retained LMO pairs to be strong.");
    }

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    if (prescreening) {
        int ijkl = 0;
        // Form unique quadruplets from the retained strong-pair connectivity graph.
        for (int ij = 0; ij < n_lmo_pairs; ij++) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            if (i > j) continue;
            for (int k : lmopair_to_lmos_[ij]) {
                if (i > k || j > k) continue;
                if (i == j && j == k) continue;

                for (int l : lmopair_to_lmos_[ij]) {
                    int kl = i_j_to_ij_[k][l];
                    if (kl == -1) continue;
                    if (i > l || j > l || k > l) continue;
                    if (i == j && j == l || j == k && k == l || i == k && k == l) continue;

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
        std::unordered_map<size_t, int> i_j_k_l_to_ijkl_new;
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
                const size_t ijkl_dense = quadruplet_key(i, j, k, l, naocc);
                if (i_j_k_l_to_ijkl_.count(ijkl_dense)) {
                    int ijkl = i_j_k_l_to_ijkl_[ijkl_dense];
                    ijkl_to_i_j_k_l_full_.push_back(std::make_tuple(i, j, k, l));
                }
            } // end l
        } // end k
    } // ij

    timer_off("Quadruples Sparsity");
}

void DLPNOCCSDT_Q::sort_quadruplets() {
    timer_on("Sort Quadruplets");

    outfile->Printf("  ==> Sorting Quadruplets <== \n\n");

    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    std::vector<std::pair<int, double>> ijkl_e_pairs(n_lmo_quadruplets);

    if (n_lmo_quadruplets == 0) {
        is_strong_quadruplet_.clear();
        qno_cutoff_.clear();
        outfile->Printf("    Number of Strong Quadruplets: %6d, Total Quadruplets: %6d, Ratio: %.4f\n\n",
                        0, 0, 0.0);
        timer_off("Sort Quadruplets");
        return;
    }

#pragma omp parallel for
    for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        ijkl_e_pairs[ijkl] = std::make_pair(ijkl, e_ijkl_[ijkl]);
    }

    std::sort(ijkl_e_pairs.begin(), ijkl_e_pairs.end(), [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return (std::fabs(a.second) > std::fabs(b.second));
    });

    double e_magnitude_total = 0.0;
    for (const auto& ijkl_e : ijkl_e_pairs) {
        e_magnitude_total += std::fabs(ijkl_e.second);
    }

    double e_magnitude_curr = 0.0;
    const bool full_quadruples_follow = algorithm_ == DLPNOMethod::CCSDTQ;
    const double t_cut_qno_full = options_.get_double("T_CUT_QNO_FULL");
    const double strong_cutoff =
        full_quadruples_follow ? t_cut_qno_full : options_.get_double("T_CUT_QNO_STRONG");
    const double weak_cutoff =
        full_quadruples_follow ? t_cut_qno_full : options_.get_double("T_CUT_QNO_WEAK");
    is_strong_quadruplet_.resize(n_lmo_quadruplets, false);
    qno_cutoff_.assign(n_lmo_quadruplets, weak_cutoff);

    int strong_count = 0;
    for (int idx = 0; idx < n_lmo_quadruplets; ++idx) {
        is_strong_quadruplet_[ijkl_e_pairs[idx].first] = true;
        qno_cutoff_[ijkl_e_pairs[idx].first] = strong_cutoff;
        e_magnitude_curr += std::fabs(ijkl_e_pairs[idx].second);
        ++strong_count;
        // For an exactly zero energy distribution, retain the previous
        // conservative behavior and classify every quadruplet as strong.
        if (e_magnitude_total > 0.0 && e_magnitude_curr / e_magnitude_total >= 0.9) break;
    }

    outfile->Printf("    Number of Strong Quadruplets: %6d, Total Quadruplets: %6d, Ratio: %.4f\n\n", strong_count, n_lmo_quadruplets, 
                            (double) strong_count / n_lmo_quadruplets);

    timer_off("Sort Quadruplets");
}

void DLPNOCCSDT_Q::qno_transform(double t_cut_qno, bool use_tuple_cutoffs) {
    timer_on("QNO transform");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    int min_qnos = options_.get_int("MIN_QNOS");
    const size_t nocc = static_cast<size_t>(std::max(naocc, 0));
    // Allowed unique occupied-index patterns are (1+1+1+1), (2+2), and
    // (2+1+1): C(N,4) + C(N,2) + 3 C(N,3). Triply repeated indices are excluded.
    const size_t n_total_possible = (nocc >= 4 ? nocc * (nocc - 1) * (nocc - 2) * (nocc - 3) / 24 : 0) +
                                    (nocc >= 2 ? nocc * (nocc - 1) / 2 : 0) +
                                    (nocc >= 3 ? nocc * (nocc - 1) * (nocc - 2) / 2 : 0);

    if (use_tuple_cutoffs && qno_cutoff_.size() != n_lmo_quadruplets) {
        throw PSIEXCEPTION(
            "Quadruplet-specific QNO cutoffs were not initialized before the iterative (Q) transform.");
    }

    double t_cut_qno_diag_scale = options_.get_double("T_CUT_QNO_DIAG_SCALE");

    X_qno_.clear();
    e_qno_.clear();
    n_qno_.clear();

    X_qno_.resize(n_lmo_quadruplets);
    e_qno_.resize(n_lmo_quadruplets);
    n_qno_.resize(n_lmo_quadruplets);

    ijkl_scale_.resize(n_lmo_quadruplets, 1.0);

    if (n_lmo_quadruplets == 0) {
        sorted_quadruplets_.clear();
        outfile->Printf("  \n");
        outfile->Printf("    Number of (Unique) Local MO quadruplets: 0\n");
        outfile->Printf("    Max Number of Possible (Unique) LMO Quadruplets: %zu (Ratio: %.4f)\n",
                        n_total_possible, 0.0);
        outfile->Printf("    No surviving quadruplets; no QNO transformation is required.\n  \n");
        timer_off("QNO transform");
        return;
    }

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

        // The pair-density sum is positive semidefinite. Treat negative
        // eigenvalues as numerical noise for selection and trace statistics;
        // using their absolute values would retain the wrong leading columns.
        double occ_total = 0.0;
        for (size_t a = 0; a < nvir_ijkl; ++a) {
            occ_total += std::max(0.0, qno_occ.get(a));
        }

        double qno_cutoff = use_tuple_cutoffs ? qno_cutoff_[ijkl] : t_cut_qno;
        if (i == j || i == k || i == l || j == k || j == l || k == l) qno_cutoff *= t_cut_qno_diag_scale;

        const size_t minimum_qnos =
            std::min(nvir_ijkl, static_cast<size_t>(std::max(1, min_qnos)));
        size_t nvir_ijkl_final = minimum_qnos;
        while (nvir_ijkl_final < nvir_ijkl && qno_occ.get(nvir_ijkl_final) >= qno_cutoff) {
            ++nvir_ijkl_final;
        }

        double occ_curr = 0.0;
        for (size_t a = 0; a < nvir_ijkl_final; ++a) {
            occ_curr += std::max(0.0, qno_occ.get(a));
        }

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(static_cast<int>(nvir_ijkl_final));

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
        occ_qno[ijkl] = std::max(0.0, qno_occ.get(n_qno_[ijkl] - 1));
        trace_qno[ijkl] = occ_total > 0.0 ? occ_curr / occ_total : 1.0;
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

    outfile->Printf("  \n");
    outfile->Printf("    Number of (Unique) Local MO quadruplets: %d\n", n_lmo_quadruplets);
    outfile->Printf("    Max Number of Possible (Unique) LMO Quadruplets: %zu (Ratio: %.4f)\n", n_total_possible,
                    n_total_possible > 0 ? (double)n_lmo_quadruplets / n_total_possible : 0.0);
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
    size_t jacobi_next_bank_memory = 0;
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
            4 * nqno3 + 16 * nqno * nlmo + 20 * nqno2;
        energy_intermediate_memory += energy_intermediates;
        quadruples_tensor_memory += 2 * nqno4;  // Gamma and T4
        jacobi_next_bank_memory += nqno4;

        // Conservative peak for one Algorithm 2 task.  The four Gamma families
        // are now constructed and released sequentially, so use their maximum
        // lexical lifetime rather than their sum.
        const size_t df_workspace =
            8 * naux * nlmo + 4 * naux * (npao + nqno) +
            2 * naux * nlmo * nqno + 2 * naux * nqno2 + 2 * naux * naux +
            npao * npao;
        const size_t gamma_terms_12_workspace =
            10 * nlmo * nqno3 + 2 * nqno4;
        const size_t gamma_term_3_workspace =
            10 * nlmo2 + 5 * nlmo * nqno2 + 3 * nqno4;
        const size_t gamma_terms_45_workspace =
            15 * nlmo * nqno2 + 10 * nqno2 + 3 * nqno4;
        const size_t gamma_term_6_workspace =
            12 * naux * nqno2 + 10 * nqno2 + 2 * nqno4;
        const size_t gamma_contraction_workspace =
            std::max({gamma_terms_12_workspace, gamma_term_3_workspace,
                      gamma_terms_45_workspace, gamma_term_6_workspace});
        const size_t energy_contraction_workspace =
            2 * nlmo * nqno3 + 5 * nqno4 + 8 * nlmo * nqno + 8 * nqno3 +
            4 * nqno2;
        gamma_workspace_per_thread =
            std::max(gamma_workspace_per_thread,
                     df_workspace + energy_intermediates +
                         std::max(gamma_contraction_workspace,
                                  energy_contraction_workspace));

        // Algorithm 5 holds a mutable T4 copy and a residual in the target QNO
        // space.  A coupled amplitude in the largest source QNO space remains
        // live while matmul_4d_permuted() ping-pongs two mixed-dimension rank-4
        // buffers.  This is the actual lexical lifetime of the refactored four
        // GEMMs, rather than a fixed multiple of max(nqno^4).
        const size_t old_dim = max_nqno;
        const size_t new_dim = nqno;
        const size_t old2 = old_dim * old_dim;
        const size_t old3 = old2 * old_dim;
        const size_t old4 = old2 * old2;
        const size_t new2 = new_dim * new_dim;
        const size_t new3 = new2 * new_dim;
        const size_t new4 = new2 * new2;
        const size_t transform_buffers = std::max(
            {2 * old3 * new_dim,
             old3 * new_dim + old2 * new2,
             2 * old2 * new2,
             old2 * new2 + old_dim * new3,
             2 * old_dim * new3,
             old_dim * new3 + new4,
             2 * new4});
        const size_t transformed_neighbor_workspace =
            old4 + old_dim * new_dim + transform_buffers;
        const size_t residual_workspace =
            2 * nqno4 + transformed_neighbor_workspace;
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

        // T1 is retained through print_results() for the diagnostic. K_iajb_,
        // by contrast, is released before this estimate and is not counted.
        for (size_t i = 0; i < i_j_to_ij_.size(); ++i) {
            const int ii = i_j_to_ij_[i][i];
            retained_lower_rank_memory += n_pno_[ii];
        }
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
        const size_t in_core_next_generation =
            write_quadruples_intermediates_ ? 0 : jacobi_next_bank_memory;
        const size_t iteration_peak = resident_memory + in_core_next_generation +
                                      nthreads * iteration_workspace_per_thread;
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
        print_memory_line("In-core Jacobi next-T4 generation",
                          write_quadruples_intermediates_ ? 0 : jacobi_next_bank_memory);
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
                4 * nqno3 + 16 * nqno * nlmo + 20 * nqno2;
            const size_t energy_contraction_workspace =
                2 * nlmo * nqno3 + 5 * nqno4 + 8 * nlmo * nqno + 8 * nqno3 +
                4 * nqno2;
            const size_t old_dim = max_nqno;
            const size_t new_dim = nqno;
            const size_t old2 = old_dim * old_dim;
            const size_t old3 = old2 * old_dim;
            const size_t old4 = old2 * old2;
            const size_t new2 = new_dim * new_dim;
            const size_t new3 = new2 * new_dim;
            const size_t new4 = new2 * new2;
            const size_t transform_buffers = std::max(
                {2 * old3 * new_dim,
                 old3 * new_dim + old2 * new2,
                 2 * old2 * new2,
                 old2 * new2 + old_dim * new3,
                 2 * old_dim * new3,
                 old_dim * new3 + new4,
                 2 * new4});
            const size_t residual_workspace =
                2 * nqno4 + old4 + old_dim * new_dim + transform_buffers;
            iteration_workspace_per_thread =
                std::max(iteration_workspace_per_thread,
                         std::max(residual_workspace,
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

        q_iv = linalg::doublet(q_iv, X_qno_[ijkl]);
        q_jv = linalg::doublet(q_jv, X_qno_[ijkl]);
        q_kv = linalg::doublet(q_kv, X_qno_[ijkl]);
        q_lv = linalg::doublet(q_lv, X_qno_[ijkl]);

        C_DGESV_shared_factorization(
            A_solve->clone(), {q_io, q_jo, q_ko, q_lo, q_iv, q_jv, q_kv, q_lv, q_ov, q_vv});

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
        std::array<Tensor<double, 4>, 10> T_mkl_list;
        auto build_gamma_terms_12 = [&]() {
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

                if (i_idx > j_idx) continue;

                // T_mkl = T_mlk with the last two virtual axes transposed, so only
                // the ten unordered positional pairs need an explicit transform.
                const int kl_position = pair_position(i_idx, j_idx);
                T_mkl_list[kl_position] =
                    Tensor<double, 4>("T_mkl", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                T_mkl_list[kl_position].zero();

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    const size_t mkl_dense = triplet_key(m, k, l, naocc);
                    if (i_j_k_to_ijk_.count(mkl_dense)) {
                        int mkl = i_j_k_to_ijk_[mkl_dense];
                        auto S_ijkl_mkl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[mkl]);
                        S_ijkl_mkl = linalg::triplet(X_qno_[ijkl], S_ijkl_mkl, X_tno_[mkl], true, false, false);
                        auto T_mkl = matmul_3d_einsums(triples_permuter_einsums(T_iajbkc_clone_[mkl], m, k, l), 
                                                        S_ijkl_mkl, n_tno_[mkl], n_qno_[ijkl]);

                        ::memcpy(&T_mkl_list[kl_position](m_ijkl, 0, 0, 0), T_mkl.data(),
                                 nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
                    } // end if
                    } // end m_ijkl
                } // end j_idx
            } // end i_idx
        };

        // Algorithm 1 intermediates for canonical Eq. (19), term 3
        std::array<Tensor<double, 2>, 10> K_minj_list;
        std::array<Tensor<double, 3>, 4> T_mkac_list;
        auto build_gamma_term_3 = [&]() {
            for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
                for (int j_idx = i_idx; j_idx < n_occupied_positions; ++j_idx) {
                    const int ij_position = pair_position(i_idx, j_idx);
                    K_minj_list[ij_position] =
                        Tensor<double, 2>("(mi|nj)", nlmo_ijkl, nlmo_ijkl);
                    einsum(0.0, Indices{index::m, index::n},
                           &K_minj_list[ij_position], 1.0,
                           Indices{index::Q, index::m}, q_io_list[i_idx],
                           Indices{index::Q, index::n}, q_io_list[j_idx]);
                }
            }

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
        };

        // Algorithm 1 intermediates for canonical Eq. (19), term 4
        std::array<Tensor<double, 3>, 4> K_iame_list;
        std::array<Tensor<double, 2>, 10> T_ijab_list;
        // Algorithm 1 intermediates for canonical Eq. (19), term 5
        std::array<Tensor<double, 3>, 4> K_mibe_list;
        auto build_gamma_terms_45 = [&]() {
            for (int idx = 0; idx < n_occupied_positions; ++idx) {
                K_iame_list[idx] = Tensor<double, 3>(
                    "(ia|me)", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::a, index::m, index::e},
                       &K_iame_list[idx], 1.0,
                       Indices{index::Q, index::a}, q_iv_list[idx],
                       Indices{index::Q, index::m, index::e}, q_ov_ein);
            }

            for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
                int i = ijkl_list[i_idx];
                for (int j_idx = i_idx; j_idx < n_occupied_positions; ++j_idx) {
                    int j = ijkl_list[j_idx];
                    int ij = i_j_to_ij_[i][j];
                    const int ij_position = pair_position(i_idx, j_idx);

                    T_ijab_list[ij_position] = Tensor<double, 2>(
                        "T_ijab", nqno_ijkl, nqno_ijkl);

                    auto S_ijkl_ij = submatrix_rows_and_cols(
                        *S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                    S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij,
                                                X_pno_[ij], true, false, false);
                    auto T_ij = linalg::triplet(S_ijkl_ij, T_iajb_[ij],
                                                S_ijkl_ij, false, false, true);
                    ::memcpy(T_ijab_list[ij_position].data(), T_ij->get_pointer(),
                             static_cast<size_t>(nqno_ijkl) * nqno_ijkl *
                                 sizeof(double));
                }
            }

            for (int idx = 0; idx < n_occupied_positions; ++idx) {
                K_mibe_list[idx] = Tensor<double, 3>(
                    "(mi|be)", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::m, index::b, index::e},
                       &K_mibe_list[idx], 1.0,
                       Indices{index::Q, index::m}, q_io_list[idx],
                       Indices{index::Q, index::b, index::e}, q_vv_ein);
            }
        };

        // Algorithm 1 theta intermediate for canonical Eq. (19), term 6.
        // Build oriented pairs lazily and release each after its final use.  A
        // quadruplet permutation never requests a diagonal positional pair, so
        // this also avoids the four unused theta blocks in the former 16-bank.
        std::array<Tensor<double, 3>, 16> theta_Qab;
        std::array<bool, 16> theta_built{};
        auto build_theta_Qab = [&](int i_idx, int j_idx) -> Tensor<double, 3>& {
            const int ij_idx = i_idx * n_occupied_positions + j_idx;
            if (theta_built[ij_idx]) return theta_Qab[ij_idx];

            const int ij_position = pair_position(i_idx, j_idx);
            theta_Qab[ij_idx] = Tensor<double, 3>(
                "theta_Qab", naux_ijkl, nqno_ijkl, nqno_ijkl);
            theta_built[ij_idx] = true;
            if (i_idx <= j_idx) {
                einsum(0.0, Indices{index::Q, index::a, index::b},
                       &theta_Qab[ij_idx], 1.0,
                       Indices{index::Q, index::a, index::e}, q_vv_ein,
                       Indices{index::e, index::b}, T_ijab_list[ij_position]);
            } else {
                einsum(0.0, Indices{index::Q, index::a, index::b},
                       &theta_Qab[ij_idx], 1.0,
                       Indices{index::Q, index::a, index::e}, q_vv_ein,
                       Indices{index::b, index::e}, T_ijab_list[ij_position]);
            }
            return theta_Qab[ij_idx];
        };

        // Algorithm 1 intermediates reused by the [Q] and (Q) energies, Eqs. (25)-(26)
        for (int i_idx = 0; i_idx < n_occupied_positions; ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = 0; j_idx < n_occupied_positions; ++j_idx) {
                int j = ijkl_list[j_idx];
                if (i_idx > j_idx) continue;
                int ij = i_j_to_ij_[i][j], ij_idx = pair_position(i_idx, j_idx);

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

        // The six canonical-Eq. (19) families have different column-symmetry
        // stabilizers. Evaluate one representative per orbit, then apply the
        // corresponding adjoint permutation with its exact orbit multiplicity:
        // terms 1/2: (k l), term 3: (i j)(k l), term 6: (i k)(j l).
        // Terms 4/5 have a trivial stabilizer and retain all 24 representatives.
        constexpr int gamma_terms_12 = 0;
        constexpr int gamma_term_3 = 1;
        constexpr int gamma_terms_45 = 2;
        constexpr int gamma_term_6 = 3;

        auto is_gamma_representative = [&](int i_idx, int j_idx, int k_idx, int l_idx,
                                           int orbit) {
            if (orbit == gamma_terms_12) return k_idx < l_idx;
            if (orbit == gamma_term_3) return i_idx < j_idx;
            if (orbit == gamma_term_6) return i_idx < k_idx;
            return true;
        };

        auto accumulate_gamma_orbit = [&](const Tensor<double, 4>& source,
                                          Tensor<double, 4>& permutation_buffer,
                                          size_t occupied_permutation, int orbit,
                                          double multiplicity) {
            einsums::for_sequence<24UL>([&](auto target_perm_idx) {
                auto &[target_i_idx, target_j_idx, target_k_idx, target_l_idx] =
                    quadruple_permutations_[target_perm_idx];
                if (!is_gamma_representative(target_i_idx, target_j_idx, target_k_idx,
                                             target_l_idx, orbit))
                    return;
                const size_t target_permutation =
                    quadruplet_key(ijkl_list[target_i_idx], ijkl_list[target_j_idx],
                                    ijkl_list[target_k_idx], ijkl_list[target_l_idx], naocc);
                if (target_permutation != occupied_permutation) return;

                permute(Indices{index::a, index::b, index::c, index::d},
                        &permutation_buffer, std::get<target_perm_idx>(einsum_indices), source);
                permutation_buffer *= multiplicity;
                gamma_ijkl += permutation_buffer;
            });
        };

        // Terms 1 and 2: swapping the last two T3 columns leaves each canonicalized
        // contraction invariant, reducing both families from 24 to 12 evaluations.
        build_gamma_terms_12();
        std::unordered_set<size_t> computed_terms_12;
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            if (!is_gamma_representative(i_idx, j_idx, k_idx, l_idx, gamma_terms_12)) return;
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx],
                l = ijkl_list[l_idx];
            const size_t occupied_permutation = quadruplet_key(i, j, k, l, naocc);
            if (!computed_terms_12.insert(occupied_permutation).second) return;

            Tensor<double, 4> gamma_ijkl_perm("gamma_ijkl_perm_12", nqno_ijkl, nqno_ijkl,
                                               nqno_ijkl, nqno_ijkl);
            gamma_ijkl_perm.zero();

            // Canonical Eq. (19), term 1; DLPNO Algorithm 2:
            // (i'a|be) t_{j'k'l'}^{ecd}.
            // The m=j slice of the term-2 T_mkl block is exactly the required
            // T_jkl amplitude.  Column symmetry supplies the reversed (l,k)
            // orientation through labels, eliminating a second T3 projection.
            const int j_ijkl =
                std::find(lmoquadruplet_to_lmos_[ijkl].begin(),
                          lmoquadruplet_to_lmos_[ijkl].end(), j) -
                lmoquadruplet_to_lmos_[ijkl].begin();
            const int kl_position = pair_position(k_idx, l_idx);
            Tensor<double, 3> T_jkl("T_jkl contiguous", nqno_ijkl,
                                     nqno_ijkl, nqno_ijkl);
            auto T_jkl_view =
                T_mkl_list[kl_position](j_ijkl, All, All, All);
            permute(Indices{index::e, index::c, index::d}, &T_jkl,
                    Indices{index::e, index::c, index::d}, T_jkl_view);
            if (k_idx <= l_idx) {
                einsum(0.0, Indices{index::a, index::b, index::c, index::d},
                       &gamma_ijkl_perm, 1.0, Indices{index::a, index::b, index::e},
                       energy_intermediates.K_iabe[i_idx],
                       Indices{index::e, index::c, index::d}, T_jkl);
            } else {
                einsum(0.0, Indices{index::a, index::b, index::c, index::d},
                       &gamma_ijkl_perm, 1.0,
                       Indices{index::a, index::b, index::e},
                       energy_intermediates.K_iabe[i_idx],
                       Indices{index::e, index::d, index::c}, T_jkl);
            }

            // Canonical Eq. (19), term 2; DLPNO Algorithm 2:
            // -(i'a|j'm) t_{mk'l'}^{bcd}.
            const int ij_idx = i_idx * n_occupied_positions + j_idx;
            einsum(1.0, Indices{index::a, index::b, index::c, index::d},
                   &gamma_ijkl_perm, -1.0, Indices{index::a, index::m},
                   energy_intermediates.K_iajm[ij_idx],
                   Indices{index::m, index::b, index::c, index::d},
                   T_mkl_list[kl_position]);

            Tensor<double, 4> permutation_buffer("gamma_permutation_buffer", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            accumulate_gamma_orbit(gamma_ijkl_perm, permutation_buffer, occupied_permutation,
                                   gamma_terms_12, 2.0);
        });
        for (auto& tensor : T_mkl_list) {
            tensor = Tensor<double, 4>("released", 0, 0, 0, 0);
        }

        // Term 3: simultaneous exchange of the two integral/T2 branches accounts
        // for the partner generated by (i j)(k l).
        build_gamma_term_3();
        std::unordered_set<size_t> computed_term_3;
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            if (!is_gamma_representative(i_idx, j_idx, k_idx, l_idx, gamma_term_3)) return;
            const size_t occupied_permutation =
                quadruplet_key(ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                                ijkl_list[l_idx], naocc);
            if (!computed_term_3.insert(occupied_permutation).second) return;

            Tensor<double, 4> gamma_ijkl_perm("gamma_ijkl_perm_3", nqno_ijkl, nqno_ijkl,
                                               nqno_ijkl, nqno_ijkl);
            Tensor<double, 4> gamma_ijkl_buffer_a("gamma_ijkl_buffer_a", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 4> gamma_ijkl_buffer_b("gamma_ijkl_buffer_b", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            // Canonical Eq. (19), term 3; DLPNO Algorithm 2:
            // +(mi'|nj') t_{mk'}^{ac} t_{nl'}^{bd}.
            Tensor<double, 3> gamma_term3("gamma_term3", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            const int ij_position = pair_position(i_idx, j_idx);
            einsum(0.0, Indices{index::n, index::a, index::c}, &gamma_term3, 1.0,
                   Indices{index::m, index::n}, K_minj_list[ij_position],
                   Indices{index::m, index::a, index::c}, T_mkac_list[k_idx]);
            einsum(0.0, Indices{index::a, index::c, index::b, index::d},
                   &gamma_ijkl_buffer_a, 1.0, Indices{index::n, index::a, index::c},
                   gamma_term3, Indices{index::n, index::b, index::d}, T_mkac_list[l_idx]);
            permute(Indices{index::a, index::b, index::c, index::d},
                    &gamma_ijkl_perm, Indices{index::a, index::c, index::b, index::d},
                    gamma_ijkl_buffer_a);
            accumulate_gamma_orbit(gamma_ijkl_perm, gamma_ijkl_buffer_b, occupied_permutation,
                                   gamma_term_3, 2.0);
        });
        for (auto& tensor : K_minj_list) {
            tensor = Tensor<double, 2>("released", 0, 0);
        }

        // Terms 4 and 5 have no nontrivial column stabilizer. The packed T2 list
        // still removes six duplicate pair transformations without materializing
        // a transposed amplitude.
        build_gamma_terms_45();
        std::unordered_set<size_t> computed_terms_45;
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            const size_t occupied_permutation =
                quadruplet_key(ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                                ijkl_list[l_idx], naocc);
            if (!computed_terms_45.insert(occupied_permutation).second) return;

            Tensor<double, 4> gamma_ijkl_perm("gamma_ijkl_perm_45", nqno_ijkl, nqno_ijkl,
                                               nqno_ijkl, nqno_ijkl);
            gamma_ijkl_perm.zero();
            Tensor<double, 4> gamma_ijkl_buffer_a("gamma_ijkl_buffer_a", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 4> gamma_ijkl_buffer_b("gamma_ijkl_buffer_b", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            const int kj_position = pair_position(k_idx, j_idx);

            // Canonical Eq. (19), term 4; DLPNO Algorithm 2:
            // -2(i'a|me) t_{k'j'}^{eb} t_{ml'}^{cd}.
            // Form the natural GEMM output (a,m,b), then transpose once so the
            // following contraction sees contiguous (a,b) and m dimensions.
            Tensor<double, 3> gamma_term4("gamma_term4", nqno_ijkl, nlmo_ijkl,
                                           nqno_ijkl);
            if (k_idx <= j_idx) {
                einsum(0.0, Indices{index::a, index::m, index::b}, &gamma_term4, 1.0,
                       Indices{index::a, index::m, index::e}, K_iame_list[i_idx],
                       Indices{index::e, index::b}, T_ijab_list[kj_position]);
            } else {
                einsum(0.0, Indices{index::a, index::m, index::b}, &gamma_term4, 1.0,
                       Indices{index::a, index::m, index::e}, K_iame_list[i_idx],
                       Indices{index::b, index::e}, T_ijab_list[kj_position]);
            }
            Tensor<double, 3> gamma_term4_transposed(
                "gamma_term4_transposed", nqno_ijkl, nqno_ijkl, nlmo_ijkl);
            permute(Indices{index::a, index::b, index::m},
                    &gamma_term4_transposed,
                    Indices{index::a, index::m, index::b}, gamma_term4);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d},
                   &gamma_ijkl_perm, -2.0, Indices{index::a, index::b, index::m},
                   gamma_term4_transposed,
                   Indices{index::m, index::c, index::d},
                   T_mkac_list[l_idx]);

            // Canonical Eq. (19), term 5; DLPNO Algorithm 2:
            // -2(be|mi') t_{k'j'}^{ce} t_{ml'}^{ad}.
            Tensor<double, 3> gamma_term5("gamma_term5", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            if (k_idx <= j_idx) {
                einsum(0.0, Indices{index::m, index::b, index::c}, &gamma_term5, 1.0,
                       Indices{index::m, index::b, index::e}, K_mibe_list[i_idx],
                       Indices{index::c, index::e}, T_ijab_list[kj_position]);
            } else {
                einsum(0.0, Indices{index::m, index::b, index::c}, &gamma_term5, 1.0,
                       Indices{index::m, index::b, index::e}, K_mibe_list[i_idx],
                       Indices{index::e, index::c}, T_ijab_list[kj_position]);
            }
            einsum(0.0, Indices{index::a, index::d, index::b, index::c},
                   &gamma_ijkl_buffer_a, 1.0, Indices{index::m, index::a, index::d},
                   T_mkac_list[l_idx], Indices{index::m, index::b, index::c}, gamma_term5);
            permute(Indices{index::a, index::b, index::c, index::d},
                    &gamma_ijkl_buffer_b, Indices{index::a, index::d, index::b, index::c},
                    gamma_ijkl_buffer_a);
            gamma_ijkl_buffer_b *= 2.0;
            gamma_ijkl_perm -= gamma_ijkl_buffer_b;

            accumulate_gamma_orbit(gamma_ijkl_perm, gamma_ijkl_buffer_a, occupied_permutation,
                                   gamma_terms_45, 1.0);
        });
        for (auto& tensor : K_iame_list) {
            tensor = Tensor<double, 3>("released", 0, 0, 0);
        }
        for (auto& tensor : K_mibe_list) {
            tensor = Tensor<double, 3>("released", 0, 0, 0);
        }
        for (auto& tensor : T_mkac_list) {
            tensor = Tensor<double, 3>("released", 0, 0, 0);
        }

        // Term 6: exchanging the two T2/DF branches, (i j) <-> (k l), leaves
        // the canonicalized contraction invariant.
        std::array<int, 16> theta_use_count{};
        std::unordered_set<size_t> counted_term_6;
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            if (!is_gamma_representative(i_idx, j_idx, k_idx, l_idx, gamma_term_6)) return;
            const size_t occupied_permutation =
                quadruplet_key(ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                                ijkl_list[l_idx], naocc);
            if (!counted_term_6.insert(occupied_permutation).second) return;
            ++theta_use_count[i_idx * n_occupied_positions + j_idx];
            ++theta_use_count[k_idx * n_occupied_positions + l_idx];
        });
        std::unordered_set<size_t> computed_term_6;
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            if (!is_gamma_representative(i_idx, j_idx, k_idx, l_idx, gamma_term_6)) return;
            const size_t occupied_permutation =
                quadruplet_key(ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                                ijkl_list[l_idx], naocc);
            if (!computed_term_6.insert(occupied_permutation).second) return;

            Tensor<double, 4> gamma_ijkl_perm("gamma_ijkl_perm_6", nqno_ijkl, nqno_ijkl,
                                               nqno_ijkl, nqno_ijkl);
            const int ij_idx = i_idx * n_occupied_positions + j_idx;
            const int kl_idx = k_idx * n_occupied_positions + l_idx;
            // Canonical Eq. (19), term 6; DLPNO Algorithm 2:
            // +(cf|ae) t_{i'j'}^{eb} t_{k'l'}^{fd}.
            auto& theta_ij = build_theta_Qab(i_idx, j_idx);
            auto& theta_kl = build_theta_Qab(k_idx, l_idx);
            einsum(0.0, Indices{index::a, index::b, index::c, index::d},
                   &gamma_ijkl_perm, 1.0, Indices{index::Q, index::a, index::b},
                   theta_ij, Indices{index::Q, index::c, index::d}, theta_kl);
            Tensor<double, 4> permutation_buffer("gamma_permutation_buffer", nqno_ijkl,
                                                  nqno_ijkl, nqno_ijkl, nqno_ijkl);
            accumulate_gamma_orbit(gamma_ijkl_perm, permutation_buffer, occupied_permutation,
                                   gamma_term_6, 2.0);

            if (--theta_use_count[ij_idx] == 0) {
                theta_Qab[ij_idx] = Tensor<double, 3>("released", 0, 0, 0);
                theta_built[ij_idx] = false;
            }
            if (--theta_use_count[kl_idx] == 0) {
                theta_Qab[kl_idx] = Tensor<double, 3>("released", 0, 0, 0);
                theta_built[kl_idx] = false;
            }
        });
        for (auto& tensor : T_ijab_list) {
            tensor = Tensor<double, 2>("released", 0, 0);
        }
        for (auto& tensor : theta_Qab) {
            tensor = Tensor<double, 3>("released", 0, 0, 0);
        }

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
    std::unordered_map<size_t, double> e_perm_energy;

    // Materialize an oriented view of a packed pair tensor only when the requested
    // occupied order is reversed. Column symmetry supplies X_ji(a,b) = X_ij(b,a).
    auto oriented_pair_tensor = [&](const auto& packed, int first_idx, int second_idx,
                                    const std::string& name,
                                    Tensor<double, 2>& transpose_workspace)
        -> const Tensor<double, 2>& {
        const Tensor<double, 2>& stored = packed[pair_position(first_idx, second_idx)];
        if (first_idx <= second_idx) return stored;
        transpose_workspace = Tensor<double, 2>(name, nqno_ijkl, nqno_ijkl);
        permute(Indices{index::a, index::b}, &transpose_workspace,
                Indices{index::b, index::a}, stored);
        return transpose_workspace;
    };

    // Assign each distinct occupied permutation to the first positional-pair
    // representative encountered in the historical 24-permutation order.
    // Processing that schedule by pair lets one T_ijm block serve all of its
    // permutations and then be released; repeated occupied indices retain the
    // exact same representative as before.
    std::array<std::vector<size_t>, 10> pair_permutation_schedule;
    std::unordered_set<size_t> scheduled_permutations;
    for (size_t perm_idx = 0; perm_idx < quadruple_permutations_.size(); ++perm_idx) {
        const auto& [i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
        const size_t occupied_permutation =
            quadruplet_key(ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                            ijkl_list[l_idx], naocc);
        if (!scheduled_permutations.insert(occupied_permutation).second) continue;
        pair_permutation_schedule[pair_position(i_idx, j_idx)].push_back(perm_idx);
    }

    for (int scheduled_pair = 0; scheduled_pair < 10; ++scheduled_pair) {
        if (pair_permutation_schedule[scheduled_pair].empty()) continue;

        const auto& first_permutation =
            quadruple_permutations_[pair_permutation_schedule[scheduled_pair].front()];
        const int pair_i_idx = std::min(std::get<0>(first_permutation),
                                        std::get<1>(first_permutation));
        const int pair_j_idx = std::max(std::get<0>(first_permutation),
                                        std::get<1>(first_permutation));
        const int pair_i = ijkl_list[pair_i_idx];
        const int pair_j = ijkl_list[pair_j_idx];

        Tensor<double, 4> T_ijm_pair("T_ijm", nlmo_ijkl, nqno_ijkl,
                                      nqno_ijkl, nqno_ijkl);
        T_ijm_pair.zero();
        for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
            const int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
            const size_t ijm_dense = triplet_key(pair_i, pair_j, m, naocc);
            if (!i_j_k_to_ijk_.count(ijm_dense)) continue;
            const int ijm = i_j_k_to_ijk_[ijm_dense];
            auto S_ijkl_ijm = submatrix_rows_and_cols(
                *S_pao_, lmoquadruplet_to_paos_[ijkl], lmotriplet_to_paos_[ijm]);
            S_ijkl_ijm = linalg::triplet(X_qno_[ijkl], S_ijkl_ijm,
                                          X_tno_[ijm], true, false, false);
            auto T_ijm = matmul_3d_einsums(
                triples_permuter_einsums(T_iajbkc_clone_[ijm], pair_i, pair_j, m),
                S_ijkl_ijm, n_tno_[ijm], nqno_ijkl);
            ::memcpy(&T_ijm_pair(m_ijkl, 0, 0, 0), T_ijm.data(),
                     static_cast<size_t>(nqno_ijkl) * nqno_ijkl * nqno_ijkl *
                         sizeof(double));
        }

        // The schedule contains both orientations of this unordered occupied
        // pair. Pack the reversed column order once; using a contiguous tensor
        // here lets every sixth-order contraction below be a GEMM rather than
        // repeatedly operating on a logically transposed rank-four view.
        Tensor<double, 4> T_jim_pair("T_jim", nlmo_ijkl, nqno_ijkl,
                                      nqno_ijkl, nqno_ijkl);
        permute(Indices{index::m, index::a, index::b, index::c}, &T_jim_pair,
                Indices{index::m, index::b, index::a, index::c}, T_ijm_pair);

        for (size_t perm_idx : pair_permutation_schedule[scheduled_pair]) {
        auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
        int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
        const size_t ijkl_idx = quadruplet_key(i, j, k, l, naocc);

        if (!e_perm_energy.count(ijkl_idx)) {
            // Set up e_perm_energy
            e_perm_energy[ijkl_idx] = 0.0;

            // Get quadruples amplitude
            Tensor<double, 4> T_ijkl = quadruples_permuter(T4, i, j, k, l);
            const Tensor<double, 4>& T_ijm_oriented =
                (i_idx <= j_idx) ? T_ijm_pair : T_jim_pair;

            // Fifth-order [Q] contribution: canonical Eq. (25), DLPNO Algorithm 3.
            // Column symmetry gives stabilizers of order 2, 1, and 4 for the
            // three displayed contractions, respectively. For a distinct-index
            // quadruplet this reduces 72 primitive evaluations to 12 + 24 + 6 = 42.
            Tensor<double, 2> U_kl_workspace;
            const Tensor<double, 2>& U_kl = oriented_pair_tensor(
                intermediates.U_iajb, k_idx, l_idx, "U_kl", U_kl_workspace);
            Tensor<double, 2> K_ij_workspace;
            const Tensor<double, 2>& K_ij = oriented_pair_tensor(
                intermediates.K_iajb, i_idx, j_idx, "K_ij", K_ij_workspace);
            Tensor<double, 2> L_ij = K_ij;
            L_ij *= 2.0;
            Tensor<double, 2> K_ji_workspace;
            L_ij -= oriented_pair_tensor(intermediates.K_iajb, j_idx, i_idx,
                                          "K_ji", K_ji_workspace);

            const size_t q5_cross_partner = quadruplet_key(j, i, l, k, naocc);
            const bool evaluate_q5_cross = ijkl_idx <= q5_cross_partner;
            const double q5_cross_multiplicity =
                (ijkl_idx == q5_cross_partner) ? 1.0 : 2.0;

            const std::array<size_t, 4> q5_pair_orbit = {
                ijkl_idx, quadruplet_key(j, i, k, l, naocc),
                quadruplet_key(i, j, l, k, naocc), q5_cross_partner};
            size_t q5_pair_representative = q5_pair_orbit[0];
            int q5_pair_multiplicity = 0;
            for (size_t orbit_idx = 0; orbit_idx < q5_pair_orbit.size(); ++orbit_idx) {
                q5_pair_representative =
                    std::min(q5_pair_representative, q5_pair_orbit[orbit_idx]);
                bool first_occurrence = true;
                for (size_t previous = 0; previous < orbit_idx; ++previous) {
                    if (q5_pair_orbit[previous] == q5_pair_orbit[orbit_idx]) {
                        first_occurrence = false;
                        break;
                    }
                }
                if (first_occurrence) ++q5_pair_multiplicity;
            }
            const bool evaluate_q5_pair = ijkl_idx == q5_pair_representative;

            // Evaluate every fifth-order product as GEMV + dot.  Explicitly
            // permute T4 so the two contracted virtual indices form one
            // contiguous matrix dimension; this replaces the scalar v^4 loop
            // without changing the symmetry representatives or multiplicities.
            Tensor<double, 4> T_buffer("T_buffer", nqno_ijkl, nqno_ijkl,
                                        nqno_ijkl, nqno_ijkl);
            Tensor<double, 2> q5_reduced("q5_reduced", nqno_ijkl, nqno_ijkl);

            permute(Indices{index::a, index::c, index::b, index::d}, &T_buffer,
                    Indices{index::a, index::b, index::c, index::d}, T_ijkl);
            einsum(0.0, Indices{index::a, index::c}, &q5_reduced, 1.0,
                   Indices{index::a, index::c, index::b, index::d}, T_buffer,
                   Indices{index::b, index::d}, U_kl);
            e_perm_energy[ijkl_idx] -=
                2.0 * linear_algebra::dot(q5_reduced, L_ij);

            if (evaluate_q5_cross) {
                permute(Indices{index::c, index::d, index::a, index::b},
                        &T_buffer,
                        Indices{index::a, index::b, index::c, index::d}, T_ijkl);
                einsum(0.0, Indices{index::c, index::d}, &q5_reduced, 1.0,
                       Indices{index::c, index::d, index::a, index::b}, T_buffer,
                       Indices{index::a, index::b}, U_kl);
                e_perm_energy[ijkl_idx] +=
                    q5_cross_multiplicity *
                    linear_algebra::dot(q5_reduced, K_ij);
            }

            if (evaluate_q5_pair) {
                einsum(0.0, Indices{index::a, index::b}, &q5_reduced, 1.0,
                       Indices{index::a, index::b, index::c, index::d}, T_ijkl,
                       Indices{index::c, index::d}, U_kl);
                e_perm_energy[ijkl_idx] +=
                    q5_pair_multiplicity *
                    linear_algebra::dot(q5_reduced, L_ij);
            }

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
            // Precombine the operator once, then evaluate one GEMM for alpha
            // and one for beta.  The former factorization performed two large
            // contractions for each result and retained both dc permutations.
            Tensor<double, 4> T_2minus("(2-P_cd)T", nqno_ijkl, nqno_ijkl,
                                        nqno_ijkl, nqno_ijkl);
            Tensor<double, 2> alpha_ijm_buffer("alpha_ijm_buffer", nlmo_ijkl, nqno_ijkl);
            Tensor<double, 2> beta_ijm_buffer("beta_ijm_buffer", nlmo_ijkl, nqno_ijkl);

            T_2minus = T_bar;
            T_2minus *= 2.0;
            permute(Indices{index::a, index::b, index::d, index::c}, &T_buffer,
                    Indices{index::a, index::b, index::c, index::d}, T_bar);
            T_2minus -= T_buffer;
            einsum(0.0, Indices{index::m, index::d}, &alpha_ijm_buffer, 1.0,
                   Indices{index::m, index::a, index::b, index::c},
                   T_ijm_oriented,
                   Indices{index::a, index::b, index::c, index::d}, T_2minus);

            T_2minus = T_tilde;
            T_2minus *= 2.0;
            permute(Indices{index::a, index::b, index::d, index::c}, &T_buffer,
                    Indices{index::a, index::b, index::c, index::d}, T_tilde);
            T_2minus -= T_buffer;
            einsum(0.0, Indices{index::m, index::d}, &beta_ijm_buffer, 1.0,
                   Indices{index::m, index::a, index::b, index::c},
                   T_ijm_oriented,
                   Indices{index::a, index::b, index::c, index::d}, T_2minus);

            Tensor<double, 2> K_lk_T("K_lk_T", nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::d}, &K_lk_T, Indices{index::d, index::m},
                    intermediates.K_iajm[l_idx * n_occupied_positions + k_idx]);
            Tensor<double, 2> K_kl_T("K_kl_T", nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::d}, &K_kl_T, Indices{index::d, index::m},
                    intermediates.K_iajm[k_idx * n_occupied_positions + l_idx]);

            e_perm_energy[ijkl_idx] += 2.0 * (linear_algebra::dot(alpha_ijm_buffer, K_lk_T) + linear_algebra::dot(beta_ijm_buffer, K_kl_T));

            // 2 - P_{kl} contributions
            int k_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), k) - lmoquadruplet_to_lmos_[ijkl].begin();
            Tensor<double, 3> T_ijk("T_ijk contiguous", nqno_ijkl,
                                     nqno_ijkl, nqno_ijkl);
            auto T_ijk_view =
                T_ijm_oriented(k_ijkl, All, All, All);
            permute(Indices{index::a, index::b, index::e}, &T_ijk,
                    Indices{index::a, index::b, index::e}, T_ijk_view);
            int l_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), l) - lmoquadruplet_to_lmos_[ijkl].begin();
            Tensor<double, 3> T_ijl("T_ijl contiguous", nqno_ijkl,
                                     nqno_ijkl, nqno_ijkl);
            auto T_ijl_view =
                T_ijm_oriented(l_ijkl, All, All, All);
            permute(Indices{index::a, index::b, index::e}, &T_ijl,
                    Indices{index::a, index::b, index::e}, T_ijl_view);

            Tensor<double, 3> T_ijk_bar("T_ijk_bar", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 3> T_ijl_bar("T_ijl_bar", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 3> T_ijk_tilde("T_ijk_tilde", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 3> T_ijl_tilde("T_ijl_tilde", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijk_bar, 1.0,
                   Indices{index::a, index::b, index::c, index::d}, T_bar,
                   Indices{index::a, index::b, index::e}, T_ijk);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijl_bar, 1.0,
                   Indices{index::a, index::b, index::c, index::d}, T_bar,
                   Indices{index::a, index::b, index::e}, T_ijl);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijk_tilde, 1.0,
                   Indices{index::a, index::b, index::c, index::d}, T_tilde,
                   Indices{index::a, index::b, index::e}, T_ijk);
            einsum(0.0, Indices{index::e, index::c, index::d}, &T_ijl_tilde, 1.0,
                   Indices{index::a, index::b, index::c, index::d}, T_tilde,
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
        }
    }
        
    return quadruplet_energy;
}

double DLPNOCCSDT_Q::lccsdt_q_iterations() {
    timer_on("LCCSDT(Q) Iterations");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    outfile->Printf("\n  ==> Local CCSDT(Q) <==\n\n");
    if (n_lmo_quadruplets == 0) {
        outfile->Printf("    No surviving quadruplets; the iterative (Q) correction is zero.\n\n");
        timer_off("LCCSDT(Q) Iterations");
        return 0.0;
    }

    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                         Corr. Energy    Delta E     Max R     Time (s)\n");

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0;
    bool e_converged = false, r_converged = false;

    double f_cut = options_.get_double("F_CUT_Q");
    double t_cut_iter = options_.get_double("T_CUT_ITER_Q");

    std::vector<double> e_ijkl_old(n_lmo_quadruplets, 0.0);

    // Algorithm 5 is a Jacobi iteration: every off-diagonal Fock coupling in
    // iteration n must read the same immutable T4 generation.  Keep distinct
    // record labels for the disk-backed path and a second tensor bank for the
    // in-core path.  The generations are swapped only after the OpenMP update
    // loop has completed.
    std::string current_t4_label = "T4";
    std::string next_t4_label = "T4 Jacobi Next";

    while (!(e_converged && r_converged)) {
        // RMS residual for each LMO quadruplet, used to assess convergence.
        std::vector<double> R_iajbkcld_rms(n_lmo_quadruplets, 0.0);

        std::time_t time_start = std::time(nullptr);

        std::vector<Tensor<double, 4>> next_t4_bank;
        if (!write_quadruples_intermediates_) next_t4_bank.resize(n_lmo_quadruplets);

#pragma omp parallel for schedule(dynamic, 1)
        for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
            int ijkl = sorted_quadruplets_[ijkl_sorted];
            auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];

            int nqno_ijkl = n_qno_[ijkl];
            if (std::fabs(e_ijkl_[ijkl] - e_ijkl_old[ijkl]) <
                std::fabs(e_ijkl_old[ijkl] * t_cut_iter)) {
                // A skipped quadruplet is unchanged, but it must still be
                // present in the complete next Jacobi generation.
                if (write_quadruples_intermediates_) {
                    auto T_ijkl = load_quadruples_tensor(current_t4_label, ijkl);
                    save_quadruples_tensor(T_ijkl, next_t4_label, ijkl);
                } else {
                    next_t4_bank[ijkl] = T_iajbkcld_[ijkl];
                }
                continue;
            }

            // S integrals
            std::vector<int> quadruplet_ext_domain;
            for (int m = 0; m < naocc; ++m) {
                const size_t ijkm_dense = quadruplet_key(i, j, k, m, naocc);
                const size_t ijml_dense = quadruplet_key(i, j, m, l, naocc);
                const size_t imkl_dense = quadruplet_key(i, m, k, l, naocc);
                const size_t mjkl_dense = quadruplet_key(m, j, k, l, naocc);

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
                                               ? load_quadruples_tensor(current_t4_label, ijkl)
                                               : T_iajbkcld_[ijkl];
            auto transform_coupled_t4 = [&](int coupled_ijkl,
                                             const SharedMatrix& overlap,
                                             int p, int q, int r, int s) {
                if (write_quadruples_intermediates_) {
                    auto source =
                        load_quadruples_tensor(current_t4_label, coupled_ijkl);
                    return matmul_4d_permuted(
                        source, overlap, n_qno_[coupled_ijkl], n_qno_[ijkl],
                        p, q, r, s);
                }
                // The current Jacobi generation is immutable during this
                // parallel loop, so the in-core source can be passed by
                // reference without a full rank-four copy.
                return matmul_4d_permuted(
                    T_iajbkcld_[coupled_ijkl], overlap,
                    n_qno_[coupled_ijkl], n_qno_[ijkl], p, q, r, s);
            };

            // Gamma is consumed exactly once in this task.  Move its storage
            // into the residual instead of retaining a second nqno^4 copy.
            Tensor<double, 4> R_ijkl = std::move(gamma_ijkl);

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
                const size_t ijkm_dense = quadruplet_key(i, j, k, m, naocc);
                const size_t ijml_dense = quadruplet_key(i, j, m, l, naocc);
                const size_t imkl_dense = quadruplet_key(i, m, k, l, naocc);
                const size_t mjkl_dense = quadruplet_key(m, j, k, l, naocc);

                if (l != m && i_j_k_l_to_ijkl_.count(ijkm_dense) && std::fabs((*F_lmo_)(l, m)) >= f_cut) {
                    int ijkm = i_j_k_l_to_ijkl_[ijkm_dense];
                    std::vector<int> ijkm_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijkm]);
                    auto S_ijkl_ijkm = linalg::doublet(submatrix_rows(*S_ijkl, ijkm_idx_list), X_qno_[ijkm], true, false);
                    auto T_temp = transform_coupled_t4(
                        ijkm, S_ijkl_ijkm, i, j, k, m);
                    T_temp *= (*F_lmo_)(l, m);
                    R_ijkl -= T_temp;
                }

                if (k != m && i_j_k_l_to_ijkl_.count(ijml_dense) && std::fabs((*F_lmo_)(k, m)) >= f_cut) {
                    int ijml = i_j_k_l_to_ijkl_[ijml_dense];
                    std::vector<int> ijml_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[ijml]);
                    auto S_ijkl_ijml = linalg::doublet(submatrix_rows(*S_ijkl, ijml_idx_list), X_qno_[ijml], true, false);
                    auto T_temp = transform_coupled_t4(
                        ijml, S_ijkl_ijml, i, j, m, l);
                    T_temp *= (*F_lmo_)(k, m);
                    R_ijkl -= T_temp;
                }

                if (j != m && i_j_k_l_to_ijkl_.count(imkl_dense) && std::fabs((*F_lmo_)(j, m)) >= f_cut) {
                    int imkl = i_j_k_l_to_ijkl_[imkl_dense];
                    std::vector<int> imkl_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[imkl]);
                    auto S_ijkl_imkl = linalg::doublet(submatrix_rows(*S_ijkl, imkl_idx_list), X_qno_[imkl], true, false);
                    auto T_temp = transform_coupled_t4(
                        imkl, S_ijkl_imkl, i, m, k, l);
                    T_temp *= (*F_lmo_)(j, m);
                    R_ijkl -= T_temp;
                }

                if (i != m && i_j_k_l_to_ijkl_.count(mjkl_dense) && std::fabs((*F_lmo_)(i, m)) >= f_cut) {
                    int mjkl = i_j_k_l_to_ijkl_[mjkl_dense];
                    std::vector<int> mjkl_idx_list = index_list(quadruplet_ext_domain, lmoquadruplet_to_paos_[mjkl]);
                    auto S_ijkl_mjkl = linalg::doublet(submatrix_rows(*S_ijkl, mjkl_idx_list), X_qno_[mjkl], true, false);
                    auto T_temp = transform_coupled_t4(
                        mjkl, S_ijkl_mjkl, m, j, k, l);
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
                save_quadruples_tensor(T_ijkl, next_t4_label, ijkl);
            } else {
                next_t4_bank[ijkl] = std::move(T_ijkl);
            }
            
            R_iajbkcld_rms[ijkl] = std::sqrt(linear_algebra::dot(R_ijkl, R_ijkl)) / (nqno_ijkl * nqno_ijkl);
        }

        if (write_quadruples_intermediates_) {
            std::swap(current_t4_label, next_t4_label);
        } else {
            T_iajbkcld_.swap(next_t4_bank);
            // Destroy the old generation before the energy pass so the second
            // persistent bank contributes only to the Jacobi-update peak.
            next_t4_bank.clear();
            next_t4_bank.shrink_to_fit();
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
                auto T_ijkl = load_quadruples_tensor(current_t4_label, ijkl);
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

    // The rest of the implementation uses the canonical "T4" record name.
    // If the converged disk generation occupies the alternate label, publish
    // it back under that stable name before returning.
    if (write_quadruples_intermediates_ && current_t4_label != "T4") {
#pragma omp parallel for schedule(dynamic, 1)
        for (int ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
            auto T_ijkl = load_quadruples_tensor(current_t4_label, ijkl);
            save_quadruples_tensor(T_ijkl, "T4", ijkl);
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
        K_iajb_.clear();
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
        K_ovov_tno_cache_.clear();
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

    // Full quadruples always form the iterative-(Q) amplitudes used as their
    // starting point; Q0_APPROXIMATION applies only when DLPNO-CCSDT(Q) is the
    // final requested method.
    const bool q0_only = algorithm_ == DLPNOMethod::CCSDT_Q && options_.get_bool("Q0_APPROXIMATION");
    if (!q0_only) {
        // Step 3: iterative non-semicanonical (Q), Algorithm 5 and Eq. (62).
        outfile->Printf("\n\n  ==> Computing Full Iterative (Q) <==\n\n");

        const bool full_quadruples_follow = algorithm_ == DLPNOMethod::CCSDTQ;
        const double t_cut_qno_full = options_.get_double("T_CUT_QNO_FULL");
        const double t_cut_qno_strong =
            full_quadruples_follow ? t_cut_qno_full : options_.get_double("T_CUT_QNO_STRONG");
        const double t_cut_qno_weak =
            full_quadruples_follow ? t_cut_qno_full : options_.get_double("T_CUT_QNO_WEAK");
        outfile->Printf("     T_CUT_QNO set to %6.3e for strong quadruplets \n", t_cut_qno_strong);
        outfile->Printf("     T_CUT_QNO set to %6.3e for weak quadruplets   \n\n", t_cut_qno_weak);

        // Sort quadruplets into "strong" and "weak" quadruplets
        sort_quadruplets();
        qno_transform(t_cut_qno, true);
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

    const std::string quadruples_energy_label = q0_only ? "CCSDT(Q0)" : "CCSDT(Q)";
    const std::string quadruples_correction_label = q0_only ? "(Q0)" : "(Q)";

    set_scalar_variable(quadruples_energy_label + " CORRELATION ENERGY", e_ccsdt_q_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsdt_q_corr);
    set_scalar_variable(quadruples_energy_label + " TOTAL ENERGY", e_ccsdt_q_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsdt_q_total);
    set_scalar_variable(quadruples_correction_label + " CORRECTION ENERGY",
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

    const bool q0_approximation =
        algorithm_ == DLPNOMethod::CCSDT_Q && options_.get_bool("Q0_APPROXIMATION");
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

DLPNOCCSDTQ::DLPNOCCSDTQ(SharedWavefunction ref_wfn, Options& options)
    : DLPNOCCSDT_Q(ref_wfn, options),
      disk_qno_integrals_(options.get_bool("DLPNO_CCSDTQ_DISK_INTS")),
      extrapolate_t4_(options.get_bool("EXTRAPOLATE_T4")),
      quadruples_damping_ratio_(options.get_double("QUADRUPLES_DAMPING_RATIO")) {}
DLPNOCCSDTQ::~DLPNOCCSDTQ() {}

void DLPNOCCSDTQ::print_header() {
    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                   DLPNO-CCSDTQ                \n");
    outfile->Printf("                   by Andy Jiang               \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  Reference: J. Chem. Theory Comput. 22, 2825--2845 (2026)\n");
    outfile->Printf("             DOI: 10.1021/acs.jctc.5c01910\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Full-quadruples orbital thresholds:\n");
    outfile->Printf("    T_CUT_QNO_FULL              = %6.3e \n", options_.get_double("T_CUT_QNO_FULL"));
    outfile->Printf("    T_CUT_QNO_STRONG            = %6.3e \n", options_.get_double("T_CUT_QNO_FULL"));
    outfile->Printf("    T_CUT_QNO_WEAK              = %6.3e \n", options_.get_double("T_CUT_QNO_FULL"));
    outfile->Printf("    MIN_QNOS                    = %6d   \n", options_.get_int("MIN_QNOS"));
    outfile->Printf("    T_CUT_XPNO                  = %6.3e \n", options_.get_double("T_CUT_XPNO"));
    outfile->Printf("    T_CUT_XPNO_DIAG_SCALE       = %6.3e \n", options_.get_double("T_CUT_XPNO_DIAG_SCALE"));
    outfile->Printf("    T_CUT_TRACE_XPNO            = %6.3e \n", options_.get_double("T_CUT_TRACE_XPNO"));
    outfile->Printf("    MIN_PNOS                    = %6d   \n", options_.get_int("MIN_PNOS"));
    outfile->Printf("\n  Quadruplet domains and local density fitting:\n");
    outfile->Printf("    T_CUT_DO_QUADS              = %6.3e \n", options_.get_double("T_CUT_DO_QUADS"));
    outfile->Printf("    T_CUT_MKN_QUADS             = %6.3e \n", options_.get_double("T_CUT_MKN_QUADS"));
    outfile->Printf("\n  Iterative solver settings:\n");
    outfile->Printf("    E_CONVERGENCE               = %6.3e \n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE               = %6.3e \n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("    DLPNO_MAXITER               = %6d   \n", options_.get_int("DLPNO_MAXITER"));
    outfile->Printf("    DIIS_MAX_VECS               = %6d   \n", options_.get_int("DIIS_MAX_VECS"));
    outfile->Printf("    DLPNO_QUADS_MICROITERATIONS = %6d   \n", options_.get_int("DLPNO_QUADS_MICROITERATIONS"));
    outfile->Printf("    QUADRUPLES_DAMPING_RATIO     = %6.3e \n", quadruples_damping_ratio_);
    outfile->Printf("\n  Memory controls:\n");
    outfile->Printf("    DLPNO_TOGGLE_MEMORY         = %6s   \n", toggle_memory_ ? "TRUE" : "FALSE");
    outfile->Printf("    DLPNO_CCSDTQ_DISK_INTS      = %6s   \n", disk_qno_integrals_ ? "TRUE" : "FALSE");
    outfile->Printf("    EXTRAPOLATE_T4              = %6s   \n\n", extrapolate_t4_ ? "TRUE" : "FALSE");
    outfile->Printf("  These are the requested memory controls; the memory estimator may\n");
    outfile->Printf("  enable disk-backed QNO integrals or disable T4 DIIS below.\n\n");
    outfile->Printf("  Lower-rank PNO/TNO thresholds and disk settings were reported by\n");
    outfile->Printf("  the preceding DLPNO-CCSD and DLPNO-CCSDT stages.\n\n");
}

void DLPNOCCSDTQ::estimate_memory() {
    const size_t naocc = i_j_to_ij_.size();
    const size_t n_lmo_pairs = ij_to_i_j_.size();
    const size_t n_lmo_triplets = ijk_to_i_j_k_.size();
    const size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();
    const size_t nthreads = std::max<size_t>(1, static_cast<size_t>(nthread_));

    size_t max_nqno = 0;
    size_t max_npno = 0;
    size_t max_xpno = 0;
    for (const int nqno : n_qno_) max_nqno = std::max(max_nqno, static_cast<size_t>(nqno));
    for (const int npno : n_pno_) max_npno = std::max(max_npno, static_cast<size_t>(npno));
    for (const int nxpno : n_xpno_) max_xpno = std::max(max_xpno, static_cast<size_t>(nxpno));
    const size_t max_npno2 = max_npno * max_npno;
    const size_t max_xpno2 = max_xpno * max_xpno;
    const size_t max_xpno4 = max_xpno2 * max_xpno2;

    size_t pno_basis_memory = 0;
    size_t singles_words = 0;
    size_t projected_pair_singles_memory = 0;
    size_t doubles_words = 0;
    size_t triples_words = 0;
    for (size_t i = 0; i < naocc; ++i) {
        const int ii = i_j_to_ij_[i][i];
        singles_words += static_cast<size_t>(n_pno_[ii]);
    }
    for (size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        const auto& [i, j] = ij_to_i_j_[ij];
        const size_t npao = lmopair_to_paos_[ij].size();
        const size_t npno = n_pno_[ij];
        doubles_words += npno * npno;
        projected_pair_singles_memory += lmopair_to_lmos_[ij].size() * npno;
        // Opposite pair orientations share the same transformation and energy vector.
        if (i <= j) pno_basis_memory += npao * npno + npno;
    }
    for (size_t ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        const size_t ntno = n_tno_[ijk];
        triples_words += ntno * ntno * ntno;
    }

    size_t qno_basis_memory = 0;
    size_t quadruples_amplitude_memory = 0;
    size_t xpno_basis_memory = 0;
    size_t always_in_core_df_memory = 0;
    size_t disk_eligible_df_memory = 0;
    size_t projected_quadruplet_singles_memory = 0;
    size_t quadruples_residual_memory = 0;
    size_t xpno_projected_t4_memory = 0;
    size_t delta_l_memory = 0;
    size_t delta_m_memory = 0;

    size_t integral_workspace_per_thread = 0;
    size_t delta_lm_workspace_per_thread = 0;
    size_t doubles_t4_workspace_per_thread = 0;
    size_t triples_t4_workspace_per_thread = 0;
    size_t quadruples_workspace_per_thread = 0;
    size_t spin_workspace_per_thread = 0;
    size_t xpno_transform_workspace_per_thread = 0;
    size_t xpno_stream_workspace_per_thread = 0;
    size_t disk_integral_workspace_per_thread = 0;

    // XPNO transforms are shared by the two orientations of each pair. The
    // projected T_mnkl blocks implement main-text Eqs. (83)--(85).
    for (size_t kl = 0; kl < n_lmo_pairs; ++kl) {
        const auto& [k, l] = ij_to_i_j_[kl];
        const size_t npno = n_pno_[kl];
        const size_t nlmo = lmopair_to_lmos_[kl].size();
        const size_t naux = lmopair_to_ribfs_[kl].size();
        const size_t npno2 = npno * npno;
        const size_t npno4 = npno2 * npno2;

        delta_l_memory += nlmo * nlmo * npno2;
        delta_m_memory += npno4;
        const size_t pair_df_words = naux * nlmo * npno;
        const size_t old_dim = max_nqno;
        const size_t new_dim = npno;
        const size_t old2 = old_dim * old_dim;
        const size_t old3 = old2 * old_dim;
        const size_t new2 = new_dim * new_dim;
        const size_t new3 = new2 * new_dim;
        const size_t new4 = new2 * new2;
        const size_t transform_buffers = std::max(
            {2 * old3 * new_dim,
             old3 * new_dim + old2 * new2,
             2 * old2 * new2,
             old2 * new2 + old_dim * new3,
             2 * old_dim * new3,
             old_dim * new3 + new4,
             2 * new4});
        // The fused routine retains one packed DF tensor and one streamed
        // g_m^{ef} slice. Delta-L and delta-M share both before advancing n.
        const size_t delta_lm_workspace =
            pair_df_words + naux * npno + nlmo * npno2 +
            max_nqno * npno + transform_buffers + 5 * npno4;
        delta_lm_workspace_per_thread =
            std::max(delta_lm_workspace_per_thread, delta_lm_workspace);

        if (k > l) continue;
        const size_t npao_ext = lmopair_to_paos_ext_[kl].size();
        const size_t nxpno = n_xpno_[kl];
        xpno_basis_memory += npao_ext * nxpno + nxpno;

        // Conservative peak for the matrices simultaneously scoped in
        // xpno_transform(): extended-PAO overlap/Fock canonicalization, the six
        // pair-density contributions, and their PAO-to-PNO transformations.
        const size_t xpno_transform_workspace =
            10 * npao_ext * npao_ext + 8 * npao_ext * max_npno +
            16 * max_npno2;
        xpno_transform_workspace_per_thread =
            std::max(xpno_transform_workspace_per_thread, xpno_transform_workspace);

        const size_t nxpno2 = nxpno * nxpno;
        const size_t nxpno4 = nxpno2 * nxpno2;
        {
            const size_t old_dim = max_nqno;
            const size_t new_dim = nxpno;
            const size_t old2 = old_dim * old_dim;
            const size_t old3 = old2 * old_dim;
            const size_t new2 = new_dim * new_dim;
            const size_t new3 = new2 * new_dim;
            const size_t new4 = new2 * new2;
            const size_t transform_buffers = std::max(
                {2 * old3 * new_dim,
                 old3 * new_dim + old2 * new2,
                 2 * old2 * new2,
                 old2 * new2 + old_dim * new3,
                 2 * old_dim * new3,
                 old_dim * new3 + new4,
                 2 * new4});
            xpno_stream_workspace_per_thread = std::max(
                xpno_stream_workspace_per_thread,
                old_dim * new_dim + transform_buffers);
        }
        for (const int m : lmopair_to_lmos_[kl]) {
            for (const int n : lmopair_to_lmos_[kl]) {
                if (m > n) continue;
                const int mn = i_j_to_ij_[m][n];
                const size_t mnkl_index = quadruplet_key(m, n, k, l, naocc);
                if (mn == -1 || i_j_k_l_to_ijkl_.count(mnkl_index) == 0) continue;
                xpno_projected_t4_memory += nxpno4;
            }
        }
    }

    // T4-dependent triples work is separate from the inherited CCSDT residual
    // workspace because the two phases execute sequentially.
    for (size_t ijk = 0; ijk < n_lmo_triplets; ++ijk) {
        const size_t naux = lmotriplet_to_ribfs_[ijk].size();
        const size_t nlmo = lmotriplet_to_lmos_[ijk].size();
        const size_t ntno = n_tno_[ijk];
        const size_t ntno2 = ntno * ntno;
        const size_t ntno3 = ntno2 * ntno;
        const size_t t1_dressed_df =
            3 * ntno + 6 * naux * (nlmo + ntno) + 2 * naux * nlmo * ntno +
            naux * ntno2 + naux + 2 * naux * nlmo * nlmo;
        const size_t local_integrals = 3 * ntno * nlmo * nlmo + nlmo * ntno3;
        const size_t old_dim = max_nqno;
        const size_t new_dim = ntno;
        const size_t old2 = old_dim * old_dim;
        const size_t old3 = old2 * old_dim;
        const size_t new2 = new_dim * new_dim;
        const size_t new3 = new2 * new_dim;
        const size_t new4 = new2 * new2;
        const size_t transform_buffers = std::max(
            {2 * old3 * new_dim,
             old3 * new_dim + old2 * new2,
             2 * old2 * new2,
             old2 * new2 + old_dim * new3,
             2 * old_dim * new3,
             old_dim * new3 + new4,
             2 * new4});
        const size_t cross_domain_t4 =
            old_dim * new_dim + transform_buffers + 3 * new4 + 5 * ntno3;
        triples_t4_workspace_per_thread =
            std::max(triples_t4_workspace_per_thread, t1_dressed_df + local_integrals + cross_domain_t4);
    }

    for (size_t ijkl = 0; ijkl < n_lmo_quadruplets; ++ijkl) {
        const size_t naux = lmoquadruplet_to_ribfs_[ijkl].size();
        const size_t nlmo = lmoquadruplet_to_lmos_[ijkl].size();
        const size_t npao = lmoquadruplet_to_paos_[ijkl].size();
        const size_t nqno = n_qno_[ijkl];
        const size_t nlmo2 = nlmo * nlmo;
        const size_t nqno2 = nqno * nqno;
        const size_t nqno3 = nqno2 * nqno;
        const size_t nqno4 = nqno2 * nqno2;

        qno_basis_memory += npao * nqno + nqno;
        quadruples_amplitude_memory += nqno4;
        quadruples_residual_memory += nqno4;
        projected_quadruplet_singles_memory += nlmo * nqno;
        always_in_core_df_memory += 4 * naux * (nlmo + nqno);
        disk_eligible_df_memory += naux * nlmo * nqno + naux * nqno2;

        // Local Matrix objects coexist briefly with the destination tensors in
        // compute_integrals(); include both PAO- and QNO-basis transforms.
        const size_t largest_solve_rhs =
            std::max({naux * nlmo, naux * nqno, naux * nlmo * nqno,
                      naux * nqno2});
        const size_t integral_workspace =
            4 * naux * (nlmo + std::max(npao, nqno)) + naux * nlmo * nqno +
            naux * nqno2 + 2 * naux * naux + nlmo * npao + npao * npao +
            npao * nqno + nqno2 + largest_solve_rhs;
        integral_workspace_per_thread = std::max(integral_workspace_per_thread, integral_workspace);
        // With disk-backed QNO integrals, the Matrix and destination Einsums
        // tensor coexist while the record is copied and saved. In-core storage
        // already includes the destination tensors in the resident count.
        const size_t disk_staging_workspace = naux * nlmo * nqno + naux * nqno2;
        disk_integral_workspace_per_thread =
            std::max(disk_integral_workspace_per_thread,
                     integral_workspace + disk_staging_workspace);

        doubles_t4_workspace_per_thread =
            std::max(doubles_t4_workspace_per_thread,
                     3 * nqno4 + 4 * nqno2 + nqno * max_npno);
        spin_workspace_per_thread = std::max(spin_workspace_per_thread, 8 * nqno4);

        // Conservative upper bound for the simultaneously live tensors in the
        // A--M construction and residual assembly (canonical Eqs. (37)--(50);
        // local SI Eqs. (S1)--(S34)). The largest groups are named after their role
        // so additions to the contraction code can be mirrored here directly.
        const size_t dressed_df_workspace =
            4 * nqno + 8 * naux * (nlmo + nqno) + 4 * naux * nlmo * nqno +
            2 * naux * nqno2 + naux + 4 * naux * nlmo2;
        const size_t projected_amplitude_workspace =
            16 * nqno2 + 8 * nlmo * nqno2 + nlmo2 * nqno2 +
            4 * nlmo2 * nqno3 + 4 * nlmo * nqno4;
        const size_t four_index_integral_workspace =
            16 * nqno * nlmo2 + 4 * nlmo2 * nqno2 + 8 * nlmo * nqno2 +
            5 * nlmo * nqno3;
        const size_t a_through_m_workspace =
            4 * nqno3 + 16 * nlmo * nqno + 2 * nqno2 +
            34 * nlmo * nqno2 + 16 * nlmo2 + 9 * nlmo * nqno3 +
            26 * nqno * nlmo2 + 10 * nqno4;
        const size_t residual_permutation_workspace = 8 * nqno4 + max_xpno4;
        quadruples_workspace_per_thread =
            std::max(quadruples_workspace_per_thread,
                     dressed_df_workspace + projected_amplitude_workspace +
                         four_index_integral_workspace + a_through_m_workspace +
                         residual_permutation_workspace);
    }

    const size_t lower_rank_residual_memory =
        singles_words + 2 * doubles_words + triples_words;
    const size_t thread_accumulation_buffers = nthreads * (singles_words + doubles_words);
    // The legacy CCSD estimate does not explicitly include PNO transforms,
    // T1 amplitudes, or the pair-domain T1 projections retained by CCSDTQ.
    const size_t inherited_supplemental_memory =
        pno_basis_memory + singles_words + projected_pair_singles_memory;
    auto iteration_vector_words = [&]() {
        return singles_words + doubles_words + triples_words +
               (extrapolate_t4_ ? quadruples_amplitude_memory : 0);
    };
    size_t in_core_disk_eligible_df_memory = disk_qno_integrals_ ? 0 : disk_eligible_df_memory;

    struct MemoryPeaks {
        size_t common_resident;
        size_t xpno;
        size_t integral;
        size_t iteration_resident;
        size_t spin;
        size_t r4;
        size_t r3;
        size_t r2;
        size_t diis;
        size_t iteration;
    };

    auto memory_peaks = [&]() {
        const size_t common_resident =
            ccsdt_resident_memory_doubles_ + inherited_supplemental_memory +
            qno_basis_memory + quadruples_amplitude_memory + xpno_basis_memory +
            always_in_core_df_memory + in_core_disk_eligible_df_memory;
        // XPNOs are built before the full-quadruples DF integrals. Include the
        // completed XPNO storage as a conservative upper bound for parallel tasks.
        const size_t xpno_peak =
            ccsdt_resident_memory_doubles_ + inherited_supplemental_memory +
            qno_basis_memory + quadruples_amplitude_memory + xpno_basis_memory +
            nthreads * xpno_transform_workspace_per_thread;
        const size_t iteration_resident =
            common_resident + projected_quadruplet_singles_memory + lower_rank_residual_memory +
            quadruples_residual_memory + thread_accumulation_buffers;
        const size_t effective_integral_workspace =
            disk_qno_integrals_ ? disk_integral_workspace_per_thread
                                : integral_workspace_per_thread;
        const size_t integral_peak =
            common_resident + nthreads * effective_integral_workspace;
        const size_t spin_peak = iteration_resident + nthreads * spin_workspace_per_thread;
        const size_t xpno_t4_resident =
            stream_xpno_t4_ ? 0 : xpno_projected_t4_memory;
        const size_t effective_quadruples_workspace =
            quadruples_workspace_per_thread +
            (stream_xpno_t4_ ? xpno_stream_workspace_per_thread : 0);
        const size_t r4_peak =
            iteration_resident + xpno_t4_resident + delta_l_memory +
            delta_m_memory +
            nthreads * std::max(delta_lm_workspace_per_thread,
                                effective_quadruples_workspace);
        const size_t r3_peak =
            iteration_resident +
            std::max(ccsdt_iteration_workspace_doubles_, nthreads * triples_t4_workspace_per_thread);
        const size_t r2_peak =
            iteration_resident +
            std::max(ccsd_iteration_workspace_doubles_, nthreads * doubles_t4_workspace_per_thread);
        // Besides the caller's solution/error vectors, the on-disk DIIS dot
        // products can load two stored error vectors simultaneously.
        const size_t diis_peak = iteration_resident + 4 * iteration_vector_words();
        const size_t iteration_peak =
            std::max({spin_peak, r4_peak, r3_peak, r2_peak, diis_peak});
        return MemoryPeaks{common_resident, xpno_peak, integral_peak,
                           iteration_resident, spin_peak, r4_peak, r3_peak,
                           r2_peak, diis_peak, iteration_peak};
    };

    const double DOUBLES_TO_GB = 1.0e-9 * sizeof(double);
    const double BYTES_TO_GB = 1.0e-9;
    auto print_estimate = [&](const char* title) {
        const auto peaks = memory_peaks();
        auto print_memory_line = [&](const std::string& label, size_t words) {
            outfile->Printf("    %-52s : %8.3f [GB]\n", label.c_str(), words * DOUBLES_TO_GB);
        };

        outfile->Printf("\n  ==> %s <==\n\n", title);
        print_memory_line("Retained DLPNO-CCSDT state", ccsdt_resident_memory_doubles_);
        print_memory_line("Inherited PNO/T1 storage supplement", inherited_supplemental_memory);
        print_memory_line("Inherited CCSDT residual workspace", ccsdt_iteration_workspace_doubles_);
        print_memory_line("QNO transforms and orbital energies", qno_basis_memory);
        print_memory_line("T4 amplitudes", quadruples_amplitude_memory);
        print_memory_line("XPNO transforms and orbital energies", xpno_basis_memory);
        print_memory_line("Always-in-core QNO DF integrals", always_in_core_df_memory);
        print_memory_line("In-core disk-eligible QNO DF integrals", in_core_disk_eligible_df_memory);
        print_memory_line("Projected T1 and all residual tensors",
                          projected_quadruplet_singles_memory + lower_rank_residual_memory +
                              quadruples_residual_memory);
        print_memory_line("Thread-private R1/R2 accumulation buffers", thread_accumulation_buffers);
        print_memory_line("XPNO-projected T_mnkl blocks (Eqs. (83)--(85))",
                          stream_xpno_t4_ ? 0 : xpno_projected_t4_memory);
        print_memory_line("Pair-domain delta-L/delta-M intermediates", delta_l_memory + delta_m_memory);
        print_memory_line("XPNO-build workspace per thread (" + std::to_string(nthreads) + ")",
                          xpno_transform_workspace_per_thread);
        print_memory_line("Streamed XPNO-T4 workspace per thread (" +
                              std::to_string(nthreads) + ")",
                          stream_xpno_t4_ ? xpno_stream_workspace_per_thread : 0);
        print_memory_line("Integral workspace per thread (" + std::to_string(nthreads) + ")",
                          disk_qno_integrals_ ? disk_integral_workspace_per_thread
                                              : integral_workspace_per_thread);
        print_memory_line("A--M/R4 workspace per thread (" + std::to_string(nthreads) + ")",
                          std::max(
                              delta_lm_workspace_per_thread,
                              quadruples_workspace_per_thread +
                                  (stream_xpno_t4_
                                       ? xpno_stream_workspace_per_thread
                                       : 0)));
        print_memory_line("T4-to-R3 workspace per thread (" + std::to_string(nthreads) + ")",
                          triples_t4_workspace_per_thread);
        print_memory_line("T4-to-R2 workspace per thread (" + std::to_string(nthreads) + ")",
                          doubles_t4_workspace_per_thread);
        print_memory_line("Spin transform workspace per thread (" + std::to_string(nthreads) + ")",
                          spin_workspace_per_thread);
        print_memory_line("Peak flattened/on-disk DIIS working vectors",
                          4 * iteration_vector_words());
        print_memory_line("Estimated common resident memory", peaks.common_resident);
        print_memory_line("Estimated XPNO-construction peak", peaks.xpno);
        print_memory_line("Estimated integral-build peak", peaks.integral);
        print_memory_line("Estimated iteration-resident memory", peaks.iteration_resident);
        print_memory_line("Estimated spin-transform peak", peaks.spin);
        print_memory_line("Estimated R4-construction peak", peaks.r4);
        print_memory_line("Estimated R3-construction peak", peaks.r3);
        print_memory_line("Estimated R2-construction peak", peaks.r2);
        print_memory_line("Estimated DIIS peak", peaks.diis);
        print_memory_line("Estimated DLPNO-CCSDTQ peak", peaks.iteration);
        outfile->Printf("    %-52s : %8.3f [GB]\n\n", "Total memory given", memory_ * BYTES_TO_GB);
    };

    print_estimate("DLPNO-CCSDTQ Memory Requirements");

    auto peaks = memory_peaks();
    size_t required_memory = std::max({peaks.xpno, peaks.integral, peaks.iteration});
    if (toggle_memory_ && !disk_qno_integrals_ &&
        required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total required memory is more than 90%% of available memory.\n");
        outfile->Printf(
            "    Switching (Q_{ijkl}|m_{ijkl} a_{ijkl}) and "
            "(Q_{ijkl}|a_{ijkl} b_{ijkl}) to disk...\n");
        disk_qno_integrals_ = true;
        in_core_disk_eligible_df_memory = 0;
        peaks = memory_peaks();
        required_memory = std::max({peaks.xpno, peaks.integral, peaks.iteration});
        print_estimate("Updated DLPNO-CCSDTQ Memory Requirements");
    }

    if (toggle_memory_ && !stream_xpno_t4_ &&
        required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf(
            "  Total required memory remains more than 90%% of available memory.\n");
        outfile->Printf(
            "    Switching the XPNO T_mnkl bank to block-streamed projection...\n");
        stream_xpno_t4_ = true;
        peaks = memory_peaks();
        required_memory = std::max({peaks.xpno, peaks.integral, peaks.iteration});
        print_estimate("Updated DLPNO-CCSDTQ Memory Requirements");
    }

    if (toggle_memory_ && extrapolate_t4_ &&
        required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total required memory remains more than 90%% of available memory.\n");
        outfile->Printf("    Removing T4/R4 blocks from the flattened DIIS vectors...\n");
        extrapolate_t4_ = false;
        peaks = memory_peaks();
        required_memory = std::max({peaks.xpno, peaks.integral, peaks.iteration});
        print_estimate("Updated DLPNO-CCSDTQ Memory Requirements");
    }

    if (toggle_memory_ && required_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf(
            "  Total required memory remains more than 90%% of available memory after all safe toggles.\n");
        throw PSIEXCEPTION("Too little memory given for the DLPNO-CCSDTQ algorithm.");
    }

    if (disk_qno_integrals_) {
        outfile->Printf("    Writing the largest QNO-basis DF integrals to disk.\n");
    } else {
        outfile->Printf("    Keeping all QNO-basis DF integrals in core.\n");
    }
    outfile->Printf("    XPNO-projected T4 blocks are %s.\n",
                    stream_xpno_t4_ ? "streamed on demand" : "retained through R4 construction");
    outfile->Printf("    T4 amplitudes are %sincluded in DIIS extrapolation.\n\n",
                    extrapolate_t4_ ? "" : "not ");
}

void DLPNOCCSDTQ::xpno_transform(double xpno_tolerance) {
    timer_on("XPNO transform");

    int naocc = nalpha_ - nfrzc();
    int nbf = basisset_->nbf();
    int n_lmo_pairs = ij_to_i_j_.size();
    int min_pnos = options_.get_int("MIN_PNOS");
    double xpno_diag_scale = options_.get_double("T_CUT_XPNO_DIAG_SCALE");
    double xpno_occ_tolerance = options_.get_double("T_CUT_TRACE_XPNO");

    lmopair_to_paos_ext_.resize(n_lmo_pairs);
    X_xpno_.resize(n_lmo_pairs);
    e_xpno_.resize(n_lmo_pairs);
    n_xpno_.resize(n_lmo_pairs);

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

            int mk = i_j_to_ij_[m][k], ml = i_j_to_ij_[m][l], nk = i_j_to_ij_[n][k], nl = i_j_to_ij_[n][l];
            if (mk == -1 || ml == -1 || nk == -1 || nl == -1) continue;

            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[mn]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[mk]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[ml]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[nk]);
            lmopair_to_paos_ext_[kl] = merge_lists(lmopair_to_paos_ext_[kl], lmopair_to_paos_[nl]);
        }
        lmopair_to_paos_ext_[lk] = lmopair_to_paos_ext_[kl];

        //                                               //
        // ==> Canonicalize extended PAOs of pair kl <== //
        //                                               //

        auto S_pao_kl_ext = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_ext_[kl]);
        auto F_pao_kl_ext = submatrix_rows_and_cols(*F_pao_, lmopair_to_paos_ext_[kl], lmopair_to_paos_ext_[kl]);

        SharedMatrix X_pao_kl_ext;
        SharedVector e_pao_kl_ext;
        std::tie(X_pao_kl_ext, e_pao_kl_ext) = orthocanonicalizer(S_pao_kl_ext, F_pao_kl_ext);

        F_pao_kl_ext = linalg::triplet(X_pao_kl_ext, F_pao_kl_ext, X_pao_kl_ext, true, false, false);

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

            const size_t mnkl_idx = quadruplet_key(m, n, k, l, naocc);
            int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;
            if (mnkl != -1) nvir_ijkl_avg += n_qno_[mnkl];
            
            quad_count++;
        }
        D_kl_sum->scale(1.0 / (6.0 * quad_count));

        // Diagonalization of pair density
        auto X_xpno_kl = std::make_shared<Matrix>("XPNO eigenvectors", nvir_kl_ext, nvir_kl_ext);
        Vector xpno_occupations("XPNO occupations", nvir_kl_ext);
        D_kl_sum->diagonalize(*X_xpno_kl, xpno_occupations, descending);

        // The extended pair density is positive semidefinite. Select a
        // contiguous prefix of its descending, nonnegative spectrum so a
        // negative roundoff eigenvalue in the tail cannot inflate the count.
        double occ_total = 0.0;
        for (size_t a = 0; a < nvir_kl_ext; ++a) {
            occ_total += std::max(0.0, xpno_occupations.get(a));
        }

        const size_t minimum_pnos =
            std::min(nvir_kl_ext, static_cast<size_t>(std::max(1, min_pnos)));
        size_t nxpno_kl = 0;
        double occ_curr = 0.0;

        double xpno_scale = 1.0;
        if (k == l) xpno_scale = xpno_diag_scale;

        for (size_t a = 0; a < nvir_kl_ext; ++a) {
            const bool below_trace_target =
                occ_total > 0.0 && occ_curr / occ_total < xpno_occ_tolerance;
            if (a < minimum_pnos || xpno_occupations.get(a) >= xpno_scale * xpno_tolerance ||
                below_trace_target) {
                occ_curr += std::max(0.0, xpno_occupations.get(a));
                ++nxpno_kl;
            } else {
                // Occupations are descending. Once the occupation and trace
                // tests both fail, no later orbital can re-enter the prefix.
                break;
            }
        } // end a

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(static_cast<int>(nxpno_kl));

        // This transformation gives orbitals that are orthonormal but not canonical
        X_xpno_kl = X_xpno_kl->get_block({zero, X_xpno_kl->rowspi()}, {zero, dim_final});
        xpno_occupations = xpno_occupations.get_block({zero, dim_final});

        SharedMatrix xpno_canonicalizer;
        SharedVector e_xpno_kl;
        std::tie(xpno_canonicalizer, e_xpno_kl) = canonicalizer(X_xpno_kl, F_pao_kl_ext);

        // This transformation gives orbitals that are orthonormal and canonical
        X_xpno_kl = linalg::doublet(X_xpno_kl, xpno_canonicalizer, false, false);
        X_xpno_kl = linalg::doublet(X_pao_kl_ext, X_xpno_kl, false, false);

        X_xpno_[kl] = X_xpno_kl;
        e_xpno_[kl] = e_xpno_kl;
        n_xpno_[kl] = X_xpno_kl->colspi(0);
        occ_xpno[kl] = std::max(0.0, xpno_occupations.get(n_xpno_[kl] - 1));
        trace_xpno[kl] = occ_total > 0.0 ? occ_curr / occ_total : 1.0;

        // account for symmetry
        if (k < l) {
            X_xpno_[lk] = X_xpno_kl;
            e_xpno_[lk] = e_xpno_kl;
            n_xpno_[lk] = X_xpno_kl->colspi(0);
            occ_xpno[lk] = occ_xpno[kl];
            trace_xpno[lk] = trace_xpno[kl];
        } // end if (k < l)
    }

    // Print out PNO domain information
    int xpno_count_total = 0, xpno_count_min = nbf, xpno_count_max = 0;
    double occ_number_total = 0.0, occ_number_min = 2.0, occ_number_max = 0.0;
    double trace_total = 0.0, trace_min = 1.0, trace_max = 0.0;
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        xpno_count_total += n_xpno_[ij];
        xpno_count_min = std::min(xpno_count_min, n_xpno_[ij]);
        xpno_count_max = std::max(xpno_count_max, n_xpno_[ij]);
        occ_number_total += occ_xpno[ij];
        occ_number_min = std::min(occ_number_min, occ_xpno[ij]);
        occ_number_max = std::max(occ_number_max, occ_xpno[ij]);
        trace_total += trace_xpno[ij];
        trace_min = std::min(trace_min, trace_xpno[ij]);
        trace_max = std::max(trace_max, trace_xpno[ij]);
    }

    outfile->Printf("  \n");
    outfile->Printf("    Extended pair natural orbitals (XPNOs) per local MO pair:\n");
    outfile->Printf("      Avg: %3d XPNOs \n", xpno_count_total / n_lmo_pairs);
    outfile->Printf("      Min: %3d XPNOs \n", xpno_count_min);
    outfile->Printf("      Max: %3d XPNOs \n", xpno_count_max);
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

        C_DGESV_shared_factorization(
            A_solve->clone(), {q_io, q_jo, q_ko, q_lo, q_iv, q_jv, q_kv, q_lv, q_ov, q_vv});

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

        if (disk_qno_integrals_) {
#pragma omp critical
            q_ov->save(psio_.get(), PSIF_DLPNO_QIA_QNO, Matrix::SubBlocks);

#pragma omp critical
            q_vv->save(psio_.get(), PSIF_DLPNO_QAB_QNO, Matrix::ThreeIndexLowerTriangle);

            q_ov_name << "(Q_ijkl | m a) " << (ijkl);
            q_vv_name << "(Q_ijkl | a b) " << (ijkl);

            q_ov_ijkl_[ijkl] = Tensor<double, 3>(q_ov->name(), 0, 0, 0);
            q_vv_ijkl_[ijkl] = Tensor<double, 3>(q_vv->name(), 0, 0, 0);
        }
    }
}

void DLPNOCCSDTQ::load_qia_qno(int ijkl) {
    if (disk_qno_integrals_) {
        int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        int nqno_ijkl = n_qno_[ijkl];

        std::stringstream q_ov_name;
        q_ov_name << "(Q_ijkl | m a) " << (ijkl);

        auto q_ov = std::make_shared<Matrix>(q_ov_name.str(), naux_ijkl, nlmo_ijkl * nqno_ijkl);
#pragma omp critical
        q_ov->load(psio_.get(), PSIF_DLPNO_QIA_QNO, Matrix::SubBlocks);

        q_ov_ijkl_[ijkl] = Tensor<double, 3>(q_ov->name(), naux_ijkl, nlmo_ijkl, nqno_ijkl);
        ::memcpy(q_ov_ijkl_[ijkl].data(), q_ov->get_pointer(), naux_ijkl * nlmo_ijkl * nqno_ijkl * sizeof(double));
    }
}

void DLPNOCCSDTQ::load_qab_qno(int ijkl) {
    if (disk_qno_integrals_) {
        int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        int nqno_ijkl = n_qno_[ijkl];

        std::stringstream q_vv_name;
        q_vv_name << "(Q_ijkl | a b) " << (ijkl);

        auto q_vv = std::make_shared<Matrix>(q_vv_name.str(), naux_ijkl, nqno_ijkl * nqno_ijkl);
#pragma omp critical
        q_vv->load(psio_.get(), PSIF_DLPNO_QAB_QNO, Matrix::ThreeIndexLowerTriangle);

        q_vv_ijkl_[ijkl] = Tensor<double, 3>(q_vv->name(), naux_ijkl, nqno_ijkl, nqno_ijkl);
        ::memcpy(q_vv_ijkl_[ijkl].data(), q_vv->get_pointer(), naux_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
    }
}

inline Tensor<double, 4> DLPNOCCSDTQ::form_alpha_ijkl(const Tensor<double, 4>& T_ijkl) {
    // Closed-shell spin adaptation of T4; Jiang et al., main-text Eq. (30).
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

inline Tensor<double, 4> DLPNOCCSDTQ::form_beta_ijkl(const Tensor<double, 4>& alpha_ijkl) {
    // Second spin-adapted T4 combination; Jiang et al., main-text Eq. (31).
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
    // Compose the alpha and beta tensors of main-text Eqs. (30)--(31), then
    // apply the final pair permutation needed by the closed-shell T4 equations.
    // This is the nonorthogonal spin-adapted representation whose metric
    // pseudoinverse is discussed by Matthews and Stanton, JCP 142, 064108
    // (2015), DOI: 10.1063/1.4907278.
    Tensor<double, 4> alpha = form_alpha_ijkl(X);
    Tensor<double, 4> beta = form_beta_ijkl(alpha);
    Tensor<double, 4> gamma = beta;
    gamma *= 2.0;
    gamma -= quadruples_permuter(beta, 0, 1, 3, 2);

    return gamma;
}

Tensor<double, 4> DLPNOCCSDTQ::quadruples_spin_desummation(const Tensor<double, 4> &X) {
    // Minimum-norm pseudoinverse for the linearly dependent spin-adapted T4
    // manifold: Matthews and Stanton, JCP 142, 064108 (2015), Eq. (28).
    // All integer weights below share the final denominator of 480.
    int i = 0, j = 1, k = 2, l = 3;

    // 7/96 term
    Tensor<double, 4> desummed = X;
    desummed *= 35.0;

    // 1/480 terms
    Tensor<double, 4> weight_1_permutations = quadruples_permuter(X, j, i, k, l);
    weight_1_permutations += quadruples_permuter(X, k, j, i, l);
    weight_1_permutations += quadruples_permuter(X, l, j, k, i);
    weight_1_permutations += quadruples_permuter(X, i, k, j, l);
    weight_1_permutations += quadruples_permuter(X, i, l, k, j);
    weight_1_permutations += quadruples_permuter(X, i, j, l, k);

    // 11/480 terms
    Tensor<double, 4> weight_11_permutations = quadruples_permuter(X, k, i, l, j);
    weight_11_permutations += quadruples_permuter(X, l, i, j, k);
    weight_11_permutations += quadruples_permuter(X, j, l, i, k);
    weight_11_permutations += quadruples_permuter(X, l, k, i, j);
    weight_11_permutations += quadruples_permuter(X, j, k, l, i);
    weight_11_permutations += quadruples_permuter(X, k, l, j, i);
    weight_11_permutations *= 11.0;

    // 1/32 term
    Tensor<double, 4> weight_15_permutations = quadruples_permuter(X, j, i, l, k);
    weight_15_permutations += quadruples_permuter(X, k, l, i, j);
    weight_15_permutations += quadruples_permuter(X, l, k, j, i);
    weight_15_permutations *= 15.0;

    desummed += weight_1_permutations;
    desummed += weight_11_permutations;
    desummed += weight_15_permutations;
    desummed *= 1.0 / 480.0;

    return desummed;
}

SharedVector DLPNOCCSDTQ::flatten_ccsdtq_diis(
    const std::vector<SharedMatrix>& matrices, const std::vector<Tensor<double, 4>>& rank4_tensors,
    bool include_t4) const {
    size_t total_size = 0;
    for (const auto& matrix : matrices) total_size += matrix->size();
    if (include_t4) {
        for (const auto& tensor : rank4_tensors) {
            total_size += static_cast<size_t>(tensor.dim(0)) * tensor.dim(1) * tensor.dim(2) * tensor.dim(3);
        }
    }

    auto flat = std::make_shared<Vector>("flattened LCCSDTQ DIIS vector", total_size);
    double* flat_data = flat->pointer();
    size_t offset = 0;
    for (const auto& matrix : matrices) {
        const size_t size = matrix->size();
        if (size == 0) continue;
        ::memcpy(flat_data + offset, matrix->pointer()[0], size * sizeof(double));
        offset += size;
    }
    if (include_t4) {
        for (const auto& tensor : rank4_tensors) {
            const size_t size = static_cast<size_t>(tensor.dim(0)) * tensor.dim(1) * tensor.dim(2) * tensor.dim(3);
            if (size == 0) continue;
            ::memcpy(flat_data + offset, tensor.data(), size * sizeof(double));
            offset += size;
        }
    }
    return flat;
}

void DLPNOCCSDTQ::copy_ccsdtq_diis(const SharedVector& flat, std::vector<SharedMatrix>& matrices,
                                    std::vector<Tensor<double, 4>>& rank4_tensors, bool include_t4) const {
    const double* flat_data = flat->pointer();
    size_t offset = 0;
    for (auto& matrix : matrices) {
        const size_t size = matrix->size();
        if (size == 0) continue;
        ::memcpy(matrix->pointer()[0], flat_data + offset, size * sizeof(double));
        offset += size;
    }
    if (include_t4) {
        for (auto& tensor : rank4_tensors) {
            const size_t size = static_cast<size_t>(tensor.dim(0)) * tensor.dim(1) * tensor.dim(2) * tensor.dim(3);
            if (size == 0) continue;
            ::memcpy(tensor.data(), flat_data + offset, size * sizeof(double));
            offset += size;
        }
    }
}

void DLPNOCCSDTQ::add_t4_to_doubles_residual(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb,
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

        // number of quadruplet natural orbitals in quadruplet domain
        const int nqno_ijkl = n_qno_[ijkl];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        std::unordered_set<size_t> computed_perms;

        // Loop over all possible quadruplet permutations
        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            const size_t ijkl_idx = quadruplet_key(i, j, k, l, naocc);
            int kl = i_j_to_ij_[k][l];

            if (!computed_perms.count(ijkl_idx)) {
                // T4, alpha/beta, and K obey simultaneous occupied/virtual
                // column symmetry.  Swapping (i,j) leaves the output block
                // unchanged; swapping (k,l) produces the transpose that the
                // final pair-buffer flush would add back to the same canonical
                // block.  Evaluate one representative of this Klein-four orbit
                // and retain the exact number of distinct occupied keys as its
                // multiplicity (important for repeated-index quadruplets).
                std::unordered_set<size_t> orbit_keys = {
                    quadruplet_key(i, j, k, l, naocc),
                    quadruplet_key(j, i, k, l, naocc),
                    quadruplet_key(i, j, l, k, naocc),
                    quadruplet_key(j, i, l, k, naocc)};
                for (const size_t orbit_key : orbit_keys) {
                    computed_perms.emplace(orbit_key);
                }

                // Get quadruples amplitude
                Tensor<double, 4> T_ijkl = quadruples_permuter(T_iajbkcld_[ijkl], i, j, k, l);

                // Form alpha_ijkl (main-text Eq. (30)).
                Tensor<double, 4> alpha = form_alpha_ijkl(T_ijkl);
                T_ijkl = Tensor<double, 4>("released", 0, 0, 0, 0);

                // Form beta_ijkl from alpha_ijkl (main-text Eq. (31)).
                Tensor<double, 4> beta = form_beta_ijkl(alpha);
                alpha = Tensor<double, 4>("released", 0, 0, 0, 0);

                // T4 contribution to R_kl^{cd}: canonical Eq. (33), local
                // main-text Eq. (93), and the contraction-ready SI Eq. (S16).
                Tensor<double, 2> K_iajb("K_iajb", nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::a, index::b}, &K_iajb, 1.0, Indices{index::Q, index::a}, q_iv_list_[ijkl][i_idx],
                        Indices{index::Q, index::b}, q_iv_list_[ijkl][j_idx]);

                Tensor<double, 2> R_iajb_cont("R_iajb_cont", nqno_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::c, index::d}, &R_iajb_cont, 0.25, Indices{index::a, index::b, index::c, index::d}, beta,
                        Indices{index::a, index::b}, K_iajb);
                R_iajb_cont *= static_cast<double>(orbit_keys.size());

                // Form QNO overlap integral (S_ijkl_kl)
                auto S_ijkl_kl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[kl]);
                S_ijkl_kl = linalg::triplet(X_qno_[ijkl], S_ijkl_kl, X_pno_[kl], true, false, false);

                // Copy into a Psi4 Matrix for the subsequent Matrix-based basis transformation.
                auto R_iajb_psi = std::make_shared<Matrix>(nqno_ijkl, nqno_ijkl);
                ::memcpy(R_iajb_psi->get_pointer(), R_iajb_cont.data(), nqno_ijkl * nqno_ijkl * sizeof(double));
                R_iajb_buffer[thread][kl]->add(linalg::triplet(S_ijkl_kl, R_iajb_psi, S_ijkl_kl, true, false, false));

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

void DLPNOCCSDTQ::add_t4_to_triples_residual(std::vector<SharedMatrix>& R_iajbkc) {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_triplets = ijk_to_i_j_k_.size();

    // Compute residual from CCSDT
    DLPNOCCSDT::compute_R_iajbkc(R_iajbkc);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
        int ijk = sorted_triplets_[ijk_sorted];
        auto &[i, j, k] = ijk_to_i_j_k_[ijk];

        int ntno_ijk = n_tno_[ijk];
        int naux_ijk = lmotriplet_to_ribfs_[ijk].size();
        int nlmo_ijk = lmotriplet_to_lmos_[ijk].size();

        // Read integrals when disk-backed storage is enabled.
        if (disk_ints_) {
            load_qia_tno(ijk);
            load_qab_tno(ijk);
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

        // Dress the second virtual index in this orientation; this is algebraically
        // equivalent to the manuscript convention and exposes a cheaper GEMM order.
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

        // Accumulate the three short permutations separately.  Transform each
        // complementary T4 block into the target TNO basis first; spin
        // adaptation and column permutation commute with this common four-axis
        // transform.  All following contractions are therefore direct GEMVs or
        // GEMMs in the final basis, and the (m,i,j,k) projection is shared by
        // all three short permutations.
        std::array<Tensor<double, 3>, 3> R_short_permutation;
        for (auto& residual : R_short_permutation) {
            residual = Tensor<double, 3>("R_ijk short permutation", ntno_ijk,
                                          ntno_ijk, ntno_ijk);
            residual.zero();
        }

        for (int m_ijk = 0; m_ijk < nlmo_ijk; ++m_ijk) {
            const int m = lmotriplet_to_lmos_[ijk][m_ijk];
            const size_t mijk_idx = quadruplet_key(m, i, j, k, naocc);
            const int mijk = i_j_k_l_to_ijkl_.count(mijk_idx)
                                  ? i_j_k_l_to_ijkl_[mijk_idx]
                                  : -1;

            if (mijk != -1) {
                auto S_ijk_mijk = submatrix_rows_and_cols(
                    *S_pao_, lmotriplet_to_paos_[ijk],
                    lmoquadruplet_to_paos_[mijk]);
                S_ijk_mijk = linalg::triplet(X_tno_[ijk], S_ijk_mijk,
                                             X_qno_[mijk], true, false, false);
                Tensor<double, 4> T_mijk_tno = matmul_4d(
                    T_iajbkcld_[mijk], S_ijk_mijk, n_qno_[mijk], ntno_ijk);

                Tensor<double, 1> F_me("F_me", ntno_ijk);
                F_me = F_ld(m_ijk, All);
                auto K_edb = K_ledb(m_ijk, All, All, All);
                Tensor<double, 3> K_aef("K_aef", ntno_ijk, ntno_ijk,
                                         ntno_ijk);
                permute(Indices{index::a, index::f, index::e}, &K_aef,
                        Indices{index::f, index::e, index::a}, K_edb);

                for (int perm_idx = 0;
                     perm_idx < static_cast<int>(short_perms_idx.size());
                     ++perm_idx) {
                    auto &[i_idx, j_idx, k_idx] = short_perms_idx[perm_idx];
                    const int p_i = ijk_idx[i_idx];
                    const int p_j = ijk_idx[j_idx];
                    const int p_k = ijk_idx[k_idx];
                    Tensor<double, 4> T_mijk = quadruples_permuter(
                        T_mijk_tno, m, p_i, p_j, p_k);
                    Tensor<double, 4> alpha_mijk = form_alpha_ijkl(T_mijk);
                    T_mijk = Tensor<double, 4>("released", 0, 0, 0, 0);

                    // F_me contribution: canonical Eq. (35), local main-text
                    // Eq. (94), and contraction-ready SI Eq. (S17).
                    einsum(1.0, Indices{index::a, index::b, index::c},
                           &R_short_permutation[perm_idx], 1.0 / 3.0,
                           Indices{index::e, index::a, index::b, index::c},
                           alpha_mijk, Indices{index::e}, F_me);

                    // (ae|mf) contribution: canonical Eq. (35), local
                    // main-text Eq. (95), and SI Eq. (S18). K_aef and alpha
                    // expose contiguous (fe) matrix dimensions.
                    einsum(1.0, Indices{index::a, index::b, index::c},
                           &R_short_permutation[perm_idx], 1.0,
                           Indices{index::a, index::f, index::e}, K_aef,
                           Indices{index::f, index::e, index::b, index::c},
                           alpha_mijk);
                }
            }

            for (int perm_idx = 0;
                 perm_idx < static_cast<int>(short_perms_idx.size());
                 ++perm_idx) {
                auto &[i_idx, j_idx, k_idx] = short_perms_idx[perm_idx];
                const int p_j = ijk_idx[j_idx];
                const int p_k = ijk_idx[k_idx];
                for (int n_ijk = 0; n_ijk < nlmo_ijk; ++n_ijk) {
                    const int n = lmotriplet_to_lmos_[ijk][n_ijk];
                    const size_t mnjk_idx =
                        quadruplet_key(m, n, p_j, p_k, naocc);
                    if (!i_j_k_l_to_ijkl_.count(mnjk_idx)) continue;
                    const int mnjk = i_j_k_l_to_ijkl_[mnjk_idx];

                    auto S_ijk_mnjk = submatrix_rows_and_cols(
                        *S_pao_, lmotriplet_to_paos_[ijk],
                        lmoquadruplet_to_paos_[mnjk]);
                    S_ijk_mnjk = linalg::triplet(
                        X_tno_[ijk], S_ijk_mnjk, X_qno_[mnjk], true, false,
                        false);
                    Tensor<double, 4> T_mnjk = matmul_4d_permuted(
                        T_iajbkcld_[mnjk], S_ijk_mnjk, n_qno_[mnjk],
                        ntno_ijk, m, n, p_j, p_k);
                    Tensor<double, 4> alpha_mnjk = form_alpha_ijkl(T_mnjk);
                    T_mnjk = Tensor<double, 4>("released", 0, 0, 0, 0);

                    // (me|ni) contribution: canonical Eq. (35), local
                    // main-text Eq. (96), and SI Eq. (S19).  Copy the strided
                    // tensor view into a contiguous vector for one GEMV.
                    Tensor<double, 1> K_meni("K_meni", ntno_ijk);
                    K_meni = K_menj_list[i_idx](All, m_ijk, n_ijk);
                    einsum(1.0, Indices{index::a, index::b, index::c},
                           &R_short_permutation[perm_idx], -1.0,
                           Indices{index::e, index::a, index::b, index::c},
                           alpha_mnjk, Indices{index::e}, K_meni);
                }
            }
        }

        Tensor<double, 3> R_ijk_cont("R_ijk_cont", ntno_ijk, ntno_ijk,
                                      ntno_ijk);
        R_ijk_cont.zero();
        for (int perm_idx = 0;
             perm_idx < static_cast<int>(short_perms_idx.size()); ++perm_idx) {
            auto &[i_idx, j_idx, k_idx] = short_perms_idx[perm_idx];
            R_ijk_cont += triples_permuter_einsums(
                R_short_permutation[perm_idx], i_idx, j_idx, k_idx);
        }

        C_DAXPY(ntno_ijk * ntno_ijk * ntno_ijk, 1.0, R_ijk_cont.data(), 1, R_iajbkc[ijk]->get_pointer(), 1);

        if (disk_ints_) {
            q_ov_[ijk] = Tensor<double, 3>(q_ov_[ijk].name(), 0, 0, 0);
            q_vv_[ijk] = Tensor<double, 3>(q_vv_[ijk].name(), 0, 0, 0);
        }
    }
}

// Pair-domain delta L/M of main-text Eqs. (98)--(102) and SI Eqs. (S13),(S15).
std::pair<std::vector<Tensor<double, 4>>, std::vector<Tensor<double, 4>>>
DLPNOCCSDTQ::compute_delta_L_M() {
    const int naocc = i_j_to_ij_.size();
    const size_t n_lmo_pairs = ij_to_i_j_.size();

    std::vector<Tensor<double, 4>> delta_L_ijk_abm_list(n_lmo_pairs);
    std::vector<Tensor<double, 4>> delta_M_ejk_abc_list(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int npno_ij = n_pno_[ij];

        delta_L_ijk_abm_list[ij] = Tensor<double, 4>(
            "delta_L_ijk_abm", nlmo_ij, nlmo_ij, npno_ij, npno_ij);
        delta_L_ijk_abm_list[ij].zero();
        delta_M_ejk_abc_list[ij] = Tensor<double, 4>(
            "delta_M_ejk_abc", npno_ij, npno_ij, npno_ij, npno_ij);
        delta_M_ejk_abc_list[ij].zero();

        // Pack B_Q^{me} as Q x (m,e), the right operand of every streamed
        // integral GEMM.  For fixed n, B_all^T B_n produces all
        // g_m^{ef}=(me|nf) slices in nlmo*npno^2 storage; both delta L and
        // delta M consume that slice before the next n is built.
        std::vector<SharedMatrix> q_ov_ij = QIA_PNO(ij);
        Tensor<double, 3> q_qme("B_Qme", naux_ij, nlmo_ij, npno_ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            ::memcpy(&q_qme(q_ij, 0, 0), q_ov_ij[q_ij]->pointer()[0],
                     static_cast<size_t>(nlmo_ij) * npno_ij * sizeof(double));
        }
        q_ov_ij.clear();

        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            const int n = lmopair_to_lmos_[ij][n_ij];

            Tensor<double, 2> q_n("B_Qnf", naux_ij, npno_ij);
            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                ::memcpy(&q_n(q_ij, 0), &q_qme(q_ij, n_ij, 0),
                         static_cast<size_t>(npno_ij) * sizeof(double));
            }

            Tensor<double, 3> g_mef("(me|nf)", nlmo_ij, npno_ij, npno_ij);
            einsum(0.0, Indices{index::m, index::e, index::f}, &g_mef, 1.0,
                   Indices{index::Q, index::m, index::e}, q_qme,
                   Indices{index::Q, index::f}, q_n);

            // delta L_ijk^{abm} += 0.5 (me|nf) alpha_nijk^{fabe}.
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                const int k = lmopair_to_lmos_[ij][k_ij];
                const size_t nijk_idx = quadruplet_key(n, i, j, k, naocc);
                if (!i_j_k_l_to_ijkl_.count(nijk_idx)) continue;
                const int nijk = i_j_k_l_to_ijkl_[nijk_idx];

                auto S_nijk_ij = submatrix_rows_and_cols(
                    *S_pao_, lmoquadruplet_to_paos_[nijk], lmopair_to_paos_[ij]);
                S_nijk_ij = linalg::triplet(X_qno_[nijk], S_nijk_ij,
                                             X_pno_[ij], true, false, false);
                Tensor<double, 4> T_nijk = matmul_4d_permuted(
                    T_iajbkcld_[nijk], S_nijk_ij->transpose(), n_qno_[nijk],
                    npno_ij, n, i, j, k);
                Tensor<double, 4> alpha_nijk = form_alpha_ijkl(T_nijk);

                // T_efab makes (e,f) and (a,b) contiguous matrix dimensions.
                permute(Indices{index::e, index::f, index::a, index::b}, &T_nijk,
                        Indices{index::f, index::a, index::b, index::e},
                        alpha_nijk);
                auto delta_L_k =
                    delta_L_ijk_abm_list[ij](k_ij, All, All, All);
                einsum(1.0, Indices{index::m, index::a, index::b}, &delta_L_k,
                       0.5, Indices{index::m, index::e, index::f}, g_mef,
                       Indices{index::e, index::f, index::a, index::b}, T_nijk);
            }

            // delta M_eij^{abc} += -0.5 (me|nf) alpha_nmij^{fabc}.
            for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
                const int m = lmopair_to_lmos_[ij][m_ij];
                const size_t nmij_idx = quadruplet_key(n, m, i, j, naocc);
                if (!i_j_k_l_to_ijkl_.count(nmij_idx)) continue;
                const int nmij = i_j_k_l_to_ijkl_[nmij_idx];

                auto S_nmij_ij = submatrix_rows_and_cols(
                    *S_pao_, lmoquadruplet_to_paos_[nmij], lmopair_to_paos_[ij]);
                S_nmij_ij = linalg::triplet(X_qno_[nmij], S_nmij_ij,
                                             X_pno_[ij], true, false, false);
                Tensor<double, 4> T_nmij = matmul_4d_permuted(
                    T_iajbkcld_[nmij], S_nmij_ij->transpose(), n_qno_[nmij],
                    npno_ij, n, m, i, j);
                Tensor<double, 4> alpha_nmij = form_alpha_ijkl(T_nmij);
                auto g_ef = g_mef(m_ij, All, All);
                einsum(1.0, Indices{index::e, index::a, index::b, index::c},
                       &delta_M_ejk_abc_list[ij], -0.5,
                       Indices{index::e, index::f}, g_ef,
                       Indices{index::f, index::a, index::b, index::c},
                       alpha_nmij);
            }
        }
    }

    return {std::move(delta_L_ijk_abm_list),
            std::move(delta_M_ejk_abc_list)};
}

void DLPNOCCSDTQ::form_T_mnkl_xpno() {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_pairs = ij_to_i_j_.size();

    // Project each T_mnkl into the XPNO domain of kl (main-text Eq. (83)).
    // The resulting blocks are reused by the X_ijkl construction and final
    // QNO projection of main-text Eqs. (84)--(85). Stored as [kl][mn].
    T_mnkl_xpno_.resize(n_lmo_pairs);

    // Loop over all pairs
#pragma omp parallel for schedule(dynamic, 1)
    for (int kl = 0; kl < n_lmo_pairs; ++kl) {
        auto &[k, l] = ij_to_i_j_[kl];
        if (k > l) continue;

        // number of LMOs in the quadruplet domain
        const int nlmo_kl = lmopair_to_lmos_[kl].size();
        T_mnkl_xpno_[kl].resize(n_lmo_pairs);

        for (int m_kl = 0; m_kl < nlmo_kl; ++m_kl) {
            int m = lmopair_to_lmos_[kl][m_kl];
            for (int n_kl = 0; n_kl < nlmo_kl; ++n_kl) {
                int n = lmopair_to_lmos_[kl][n_kl];
                int mn = i_j_to_ij_[m][n];
                if (m > n) continue;

                const size_t mnkl_idx = quadruplet_key(m, n, k, l, naocc);
                int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;

                if (mn == -1 || mnkl == -1) continue;
                
                auto S_mnkl_kl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[mnkl], lmopair_to_paos_ext_[kl]);
                S_mnkl_kl = linalg::triplet(X_qno_[mnkl], S_mnkl_kl, X_xpno_[kl], true, false, false);

                T_mnkl_xpno_[kl][mn] = matmul_4d_permuted(
                    T_iajbkcld_[mnkl], S_mnkl_kl->transpose(), n_qno_[mnkl],
                    n_xpno_[kl], m, n, k, l);
            } // end n_kl
        } // end m_kl
    } // end kl
}

void DLPNOCCSDTQ::compute_quadruples_residual(std::vector<Tensor<double, 4>>& R_iajbkcld) {

    int naocc = i_j_to_ij_.size();
    size_t n_lmo_quadruplets = ijkl_to_i_j_k_l_.size();

    // Index orders corresponding to the 24 occupied-orbital permutations.
    auto einsum_indices = std::make_tuple(Indices{a, b, c, d}, Indices{a, b, d, c}, Indices{a, c, b, d}, Indices{a, c, d, b}, 
        Indices{a, d, b, c}, Indices{a, d, c, b}, Indices{b, a, c, d}, Indices{b, a, d, c}, Indices{b, c, a, d}, Indices{b, c, d, a}, 
        Indices{b, d, a, c}, Indices{b, d, c, a}, Indices{c, a, b, d}, Indices{c, a, d, b}, Indices{c, b, a, d}, Indices{c, b, d, a}, 
        Indices{c, d, a, b}, Indices{c, d, b, a}, Indices{d, a, b, c}, Indices{d, a, c, b}, Indices{d, b, a, c}, Indices{d, b, c, a}, 
        Indices{d, c, a, b}, Indices{d, c, b, a});

    // Form the pair-domain L/M corrections together so both consume the same
    // streamed (me|nf) slices, SI Eqs. (S13) and (S15).
    std::vector<Tensor<double, 4>> delta_L_ijk_abm_list;
    std::vector<Tensor<double, 4>> delta_M_ejk_abc_list;
    std::tie(delta_L_ijk_abm_list, delta_M_ejk_abc_list) =
        compute_delta_L_M();

// Loop over unique quadruplets
#pragma omp parallel for schedule(dynamic, 1)
    for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
        int ijkl = sorted_quadruplets_[ijkl_sorted];
        auto &[i, j, k, l] = ijkl_to_i_j_k_l_[ijkl];
        std::vector<int> ijkl_list = {i, j, k, l};
        const int FOUR = ijkl_list.size();

        // number of LMOs in the quadruplet domain
        const int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();
        // number of auxiliary functions in the quadruplet domain
        const int naux_ijkl = lmoquadruplet_to_ribfs_[ijkl].size();
        // number of quadruplet natural orbitals in quadruplet domain
        const int nqno_ijkl = n_qno_[ijkl];

        // Bookkeeping for permutations of ijk (ijkl choose ijk)
        std::unordered_map<size_t, size_t> ijk_to_ijkl_perm_idx;

        einsums::for_sequence<24UL>([&](auto perm_idx) {
            auto &[i_idx, j_idx, k_idx, l_idx] = quadruple_permutations_[perm_idx];
            int i = ijkl_list[i_idx], j = ijkl_list[j_idx], k = ijkl_list[k_idx], l = ijkl_list[l_idx];
            const size_t ijk_idx = triplet_key(i, j, k, naocc);

            ijk_to_ijkl_perm_idx[ijk_idx] = perm_idx;
        });

        // Complementary occupied triples: omitting i, j, k, or l, respectively.
        std::array<std::tuple<int, int, int>, 4> complementary_triplets = {std::make_tuple(j, k, l),
            std::make_tuple(i, k, l), std::make_tuple(j, i, l), std::make_tuple(j, k, i)};

        // (T1-DRESS) NECESSARY FOCK MATRIX ELEMENTS AND INTEGRALS

        // List of singles amplitudes projected onto ijkl space
        std::array<Tensor<double, 1>, 4> T_i_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            int i_ijkl = std::find(lmoquadruplet_to_lmos_[ijkl].begin(), lmoquadruplet_to_lmos_[ijkl].end(), i) - lmoquadruplet_to_lmos_[ijkl].begin();
            T_i_list[i_idx] = Tensor<double, 1>("T_i", nqno_ijkl);
            ::memcpy(T_i_list[i_idx].data(), &(T_m_ijkl_[ijkl])(i_ijkl, 0), nqno_ijkl * sizeof(double));
        } // end i_idx

        // Read in expensive integrals
        // (Q | m a) integrals
        load_qia_qno(ijkl);
        const Tensor<double, 3>& q_ov = q_ov_ijkl_[ijkl];
        // (Q | a b) integrals
        load_qab_qno(ijkl);
        const Tensor<double, 3>& q_vv = q_vv_ijkl_[ijkl];

        // T1-dress q_io and q_iv; checked orbital indices denote T1-dressed quantities.
        std::array<Tensor<double, 2>, 4> q_io_t1_list; // (Q | m \check{i}) = (Q | m i) + (Q | m a) t_{i}^{a}
        std::array<Tensor<double, 2>, 4> q_iv_t1_list; // (Q | \check{a} \check{i}) = (Q | a i) - t_{m}^{a} (Q | m i) + (Q | a b) t_{i}^{b} - t_{m}^{a} (Q | m b) t_{i}^{b}

        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            q_io_t1_list[i_idx] = q_io_list_[ijkl][i_idx]; // (Q | m i)
            einsum(1.0, Indices{index::Q, index::m}, &q_io_t1_list[i_idx], 1.0, Indices{index::Q, index::m, index::a}, q_ov, Indices{index::a}, T_i_list[i_idx]);

            q_iv_t1_list[i_idx] = q_iv_list_[ijkl][i_idx]; // (Q | a i)
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], -1.0, Indices{index::Q, index::m}, q_io_list_[ijkl][i_idx], Indices{index::m, index::a}, T_m_ijkl_[ijkl]);
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], 1.0, Indices{index::Q, index::a, index::b}, q_vv, Indices{index::b}, T_i_list[i_idx]);
            Tensor<double, 2> q_iv_t1_temp("q_iv_t1_temp", naux_ijkl, nlmo_ijkl);
            einsum(0.0, Indices{index::Q, index::m}, &q_iv_t1_temp, 1.0, Indices{index::Q, index::m, index::b}, q_ov, Indices{index::b}, T_i_list[i_idx]);
            einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_list[i_idx], -1.0, Indices{index::Q, index::m}, q_iv_t1_temp, Indices{index::m, index::a}, T_m_ijkl_[ijkl]);
        }

        // Store q_ov with transposed orbital indices so subsequent contractions can use GEMM kernels.
        Tensor<double, 3> q_vo("q_vo", naux_ijkl, nqno_ijkl, nlmo_ijkl);
        permute(Indices{index::Q, index::a, index::m}, &q_vo, Indices{index::Q, index::m, index::a}, q_ov);

        // T1-dress the second virtual index, opposite to the paper's convention, to reduce contraction cost.
        Tensor<double, 3> q_vv_t1 = q_vv; // (Q | a \check{b}) = (Q | a b) - (Q | a m) t_{m}^{b} 
        einsum(1.0, Indices{index::Q, index::a, index::b}, &q_vv_t1, -1.0, Indices{index::Q, index::a, index::m}, q_vo, Indices{index::m, index::b}, T_m_ijkl_[ijkl]);

        // Build the T1-dressed Fock intermediates.

        // Gamma_Q is used universally for J-like contractions
        Tensor<double, 1> gamma_Q("gamma_Q", naux_ijkl);
        einsum(0.0, Indices{index::Q}, &gamma_Q, 1.0, Indices{index::Q, index::m, index::e}, q_ov, Indices{index::m, index::e}, T_m_ijkl_[ijkl]);

        // F_me (this is scoped to ensure that the intermediate tensors are not persistent in memory)
        Tensor<double, 2> F_me("F_me", nlmo_ijkl, nqno_ijkl); {
            // J contractions
            einsum(0.0, Indices{index::m, index::e}, &F_me, 2.0, Indices{index::Q, index::m, index::e}, q_ov, Indices{index::Q}, gamma_Q);
            
            // K contractions (rc|ks)t_{k}^{c} ... (mf|ne) t_{n}^{f}
            Tensor<double, 3> F_me_K_temp("F_me_K_temp", naux_ijkl, nlmo_ijkl, nlmo_ijkl);
            einsum(0.0, Indices{index::Q, index::m, index::n}, &F_me_K_temp, 1.0, Indices{index::Q, index::m, index::f}, q_ov, Indices{index::n, index::f}, T_m_ijkl_[ijkl]);
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
            einsum(0.0, Indices{index::Q, index::f}, &F_mi_K_temp, 1.0, Indices{index::Q, index::n}, q_io_list_[ijkl][i_idx], Indices{index::n, index::f}, T_m_ijkl_[ijkl]);
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
            einsum(0.0, Indices{index::Q, index::a, index::n}, &F_ae_K_temp, 1.0, Indices{index::Q, index::a, index::f}, q_vv, Indices{index::n, index::f}, T_m_ijkl_[ijkl]);
            Tensor<double, 3> F_ae_K_temp2("F_ae_K_temp2", naux_ijkl, nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::Q, index::n, index::a}, &F_ae_K_temp2, Indices{index::Q, index::a, index::n}, F_ae_K_temp);
            einsum(1.0, Indices{index::a, index::e}, &F_ae, -1.0, Indices{index::Q, index::n, index::a}, F_ae_K_temp2, Indices{index::Q, index::n, index::e}, q_ov);

            // Add the F_me contribution to F_ae
            einsum(1.0, Indices{index::a, index::e}, &F_ae, -1.0, Indices{index::m, index::a}, T_m_ijkl_[ijkl], Indices{index::m, index::e}, F_me);
        }

        // Amplitude intermediates
        std::array<Tensor<double, 2>, 16> T_ij_list; // (Project T_ij amplitudes from PNO space of ij to QNO space of ijkl)
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            for (int j_idx = i_idx; j_idx < FOUR; ++j_idx) {
                int j = ijkl_list[j_idx];
                int ij = i_j_to_ij_[i][j];

                T_ij_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("T_ij", nqno_ijkl, nqno_ijkl);
                auto S_ijkl_ij = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij, X_pno_[ij], true, false, false);

                auto T_ij = linalg::triplet(S_ijkl_ij, T_iajb_[ij], S_ijkl_ij, false, false, true);
                ::memcpy(T_ij_list[i_idx * FOUR + j_idx].data(), T_ij->get_pointer(), nqno_ijkl * nqno_ijkl * sizeof(double));

                if (i_idx != j_idx) {
                    T_ij_list[j_idx * FOUR + i_idx] = Tensor<double, 2>(
                        "T_ji", nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::b, index::a},
                            &T_ij_list[j_idx * FOUR + i_idx],
                            Indices{index::a, index::b},
                            T_ij_list[i_idx * FOUR + j_idx]);
                }
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
                    const size_t mni_idx = triplet_key(m, n, i, naocc);
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

        // T_mni_list[j](m,i,...) is exactly T_mij.  Keep that five-index bank
        // as the single source of truth and expose the needed block as a view,
        // eliminating sixteen repeated T3 projections and the resident
        // T_mij/Z_mij tensor banks.
        std::array<int, 4> occupied_domain_positions;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            occupied_domain_positions[i_idx] =
                std::find(lmoquadruplet_to_lmos_[ijkl].begin(),
                          lmoquadruplet_to_lmos_[ijkl].end(), ijkl_list[i_idx]) -
                lmoquadruplet_to_lmos_[ijkl].begin();
        }
        auto T_mij_view = [&](int i_idx, int j_idx) {
            return T_mni_list[j_idx](All, occupied_domain_positions[i_idx], All,
                                      All, All);
        };
        auto form_Z_mij = [&](int i_idx, int j_idx) {
            auto T_mij = T_mij_view(i_idx, j_idx);
            Tensor<double, 4> Z_mij("Z_mij", nlmo_ijkl, nqno_ijkl,
                                     nqno_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::a, index::b, index::c}, &Z_mij,
                    Indices{index::m, index::a, index::b, index::c}, T_mij);
            Z_mij *= 2.0;
            Tensor<double, 4> permutation_buffer(
                "Z_mij permutation", nlmo_ijkl, nqno_ijkl, nqno_ijkl,
                nqno_ijkl);
            permute(Indices{index::m, index::a, index::b, index::c},
                    &permutation_buffer,
                    Indices{index::m, index::b, index::a, index::c}, T_mij);
            Z_mij -= permutation_buffer;
            permute(Indices{index::m, index::a, index::b, index::c},
                    &permutation_buffer,
                    Indices{index::m, index::c, index::b, index::a}, T_mij);
            Z_mij -= permutation_buffer;
            return Z_mij;
        };

        // Quadruples amplitude intermediates
        // (i -> jkl), (j -> ikl), (k -> ijl), (l -> ijk)
        std::array<Tensor<double, 5>, 4> T_nijk_complement_unsorted;

        for (int idx = 0; idx < FOUR; ++idx) {
            auto &[i, j, k] = complementary_triplets[idx];
            
            T_nijk_complement_unsorted[idx] = Tensor<double, 5>("T_nijk_unsorted", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            T_nijk_complement_unsorted[idx].zero();

            for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];
                const size_t nijk_idx = quadruplet_key(n, i, j, k, naocc);
                int nijk = i_j_k_l_to_ijkl_.count(nijk_idx) ? (i_j_k_l_to_ijkl_[nijk_idx]) : -1;
                if (nijk == -1) continue;

                auto S_ijkl_nijk = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmoquadruplet_to_paos_[nijk]);
                S_ijkl_nijk = linalg::triplet(X_qno_[ijkl], S_ijkl_nijk, X_qno_[nijk], true, false, false);

                Tensor<double, 4> T_nijk_unsorted = matmul_4d(T_iajbkcld_[nijk], S_ijkl_nijk, n_qno_[nijk], n_qno_[ijkl]);
                ::memcpy(&T_nijk_complement_unsorted[idx](n_ijkl, 0, 0, 0, 0), T_nijk_unsorted.data(), nqno_ijkl * nqno_ijkl * nqno_ijkl * nqno_ijkl * sizeof(double));
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

        // A_ej^{ab}: canonical Eq. (37), local SI Eq. (S1); dimensions 4 * (e, a, b).
        std::array<Tensor<double, 3>, 4> A_ejab_list;
        for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
            int j = ijkl_list[j_idx];

            // First term of SI Eq. (S1): (ae|bj).
            A_ejab_list[j_idx] = Tensor<double, 3>("A_ejab", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], 1.0, Indices{index::Q, index::e, index::a}, 
                    q_vv_t1, Indices{index::Q, index::b}, q_iv_t1_list[j_idx]);

            // Second term of SI Eq. (S1): (me|nj) T_mn^{ab}.
            einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], 1.0, Indices{index::e, index::m, index::n}, g_menj_list[j_idx],
                    Indices{index::m, index::n, index::a, index::b}, T_mn);

            // Third term of SI Eq. (S1): 1/2 [2(mf|ae) - (me|af)] U_mj^{fb}.
            einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], 0.5, Indices{index::m, index::f, index::e, index::a},
                    L_mfae, Indices{index::m, index::f, index::b}, U_mi_list[j_idx]);

            // Fourth term of SI Eq. (S1): -(1/2 + P_ab)(me|af) T_mj^{bf}.
            // Form its direct and virtual-index-permuted contributions.
            Tensor<double, 3> A_exchange_direct("A_exchange_direct", nqno_ijkl, nqno_ijkl, nqno_ijkl);
            Tensor<double, 3> A_exchange_permuted("A_exchange_permuted", nqno_ijkl, nqno_ijkl, nqno_ijkl);

            Tensor<double, 3> T_mi_t("T_mi_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::f, index::b}, &T_mi_t, Indices{index::m, index::b, index::f}, T_mi_list[j_idx]);

            einsum(0.0, Indices{index::e, index::a, index::b}, &A_exchange_direct, 1.0, Indices{index::m, index::f, index::e, index::a},
                    g_mfae_t, Indices{index::m, index::f, index::b}, T_mi_t);
            permute(Indices{index::e, index::b, index::a}, &A_exchange_permuted, Indices{index::e, index::a, index::b}, A_exchange_direct);
            A_exchange_direct *= 0.5;
            A_ejab_list[j_idx] -= A_exchange_direct;
            A_ejab_list[j_idx] -= A_exchange_permuted;

            // Fifth term of SI Eq. (S1): -(me|nf) Z_nmj^{fab}.
            // Z_{ijk}^{abc} = 2 T_{ijk}^{abc} - T_{ijk}^{bac} - T_{ijk}^{cba}
            // Form Z_{nmj} one triplet at a time to limit peak memory.

            for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];

                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    Tensor<double, 2> g_nm_slice = g_menf_t(n_ijkl, m_ijkl, All, All); // (f, e)

                    const size_t nmj_idx = triplet_key(n, m, j, naocc);
                    if (i_j_k_to_ijk_.count(nmj_idx)) {
                        // T_mni_list[j](n,m,...) is the same already projected
                        // T_nmj block used here; avoid a second TNO->QNO transform.
                        auto T_nmj = T_mni_list[j_idx](n_ijkl, m_ijkl, All, All,
                                                       All);
                        Tensor<double, 3> Z_nmj("Z_nmj", nqno_ijkl, nqno_ijkl,
                                                 nqno_ijkl);
                        permute(Indices{index::a, index::b, index::c}, &Z_nmj,
                                Indices{index::a, index::b, index::c}, T_nmj);
                        // Z_nmj^{abc} = 2 T_nmj^{abc} - T_nmj^{bac} - T_nmj^{cba}.
                        Tensor<double, 3> T_nmj_bac("T_nmj_bac", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                        permute(Indices{index::b, index::a, index::c}, &T_nmj_bac,
                                Indices{index::a, index::b, index::c}, T_nmj);
                        Tensor<double, 3> T_nmj_cba("T_nmj_cba", nqno_ijkl, nqno_ijkl, nqno_ijkl);
                        permute(Indices{index::c, index::b, index::a}, &T_nmj_cba,
                                Indices{index::a, index::b, index::c}, T_nmj);

                        Z_nmj *= 2.0;
                        Z_nmj -= T_nmj_bac;
                        Z_nmj -= T_nmj_cba;

                        // Accumulate the fifth contraction of SI Eq. (S1).
                        einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], -1.0, Indices{index::f, index::e},
                                g_nm_slice, Indices{index::f, index::a, index::b}, Z_nmj);
                    }
                } // end m_ijkl
            } // end n_ijkl

            // Final term of SI Eq. (S1): -F_me T_mj^{ab}.
            einsum(1.0, Indices{index::e, index::a, index::b}, &A_ejab_list[j_idx], -1.0, Indices{index::m, index::e}, F_me,
                    Indices{index::m, index::a, index::b}, T_mi_list[j_idx]);
        }

        // B_ij^{am}: canonical Eq. (38), local SI Eq. (S2); dimensions 16 * (m, a).
        std::array<Tensor<double, 2>, 16> B_ijam_list;

        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];

            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                int j = ijkl_list[j_idx];

                // First term of SI Eq. (S2): (ai|mj).
                B_ijam_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("B_ijam", nlmo_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::m}, q_io_t1_list[j_idx], 
                        Indices{index::Q, index::a}, q_iv_t1_list[i_idx]);

                // Second term of SI Eq. (S2): (mf|ae) T_ji^{fe}.
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::m, index::a, index::f, index::e},
                        g_mfae_u, Indices{index::f, index::e}, T_ij_list[j_idx * FOUR + i_idx]);

                // Third term of SI Eq. (S2): 1/2 [2(ne|mj) - (me|nj)] U_ni^{ea}.
                Tensor<double, 3> U_mi_t("U_mi_t", nqno_ijkl, nlmo_ijkl, nqno_ijkl);
                permute(Indices{index::e, index::n, index::a}, &U_mi_t, Indices{index::n, index::e, index::a}, U_mi_list[i_idx]);
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 0.5, Indices{index::e, index::n, index::m}, 
                        L_menj_list[j_idx], Indices{index::e, index::n, index::a}, U_mi_t);

                // Fourth/fifth terms of SI Eq. (S2): -(1/2 + P_ij)(me|nj) T_ni^{ae}.
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

                // Sixth term of SI Eq. (S2): F_me T_ij^{ae}.
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::m, index::e},
                        F_me, Indices{index::a, index::e}, T_ij_list[i_idx * FOUR + j_idx]);
                        
                // Final term of SI Eq. (S2): (me|nf) Z_nij^{fae}.
                // Transpose Z_nij into the contraction order used below.
                Tensor<double, 4> Z_nij = form_Z_mij(i_idx, j_idx);
                Tensor<double, 4> Z_nij_transposed("Z_nij_transposed", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::n, index::e, index::f, index::a},
                        &Z_nij_transposed,
                        Indices{index::n, index::f, index::a, index::e}, Z_nij);
                einsum(1.0, Indices{index::m, index::a}, &B_ijam_list[i_idx * FOUR + j_idx], 1.0, Indices{index::m, index::n, index::e, index::f}, g_menf_t,
                        Indices{index::n, index::e, index::f, index::a}, Z_nij_transposed);
            } // end j_idx
        } // end i_idx

        // The (m,n,e,f) orientation is needed only by the adjacent tilde-F
        // contractions. Release it before constructing E--M instead of
        // retaining a second full L_menf tensor through the complete residual.
        Tensor<double, 4> L_menf_t("L_mnef_t", nlmo_ijkl, nlmo_ijkl,
                                    nqno_ijkl, nqno_ijkl);
        permute(Indices{index::m, index::n, index::e, index::f}, &L_menf_t,
                Indices{index::m, index::e, index::n, index::f}, L_menf);

        // Tilde F_ae: canonical Eq. (39), local SI Eq. (S3); dimensions (a, e).
        Tensor<double, 2> F_ae_tilde = F_ae;
        einsum(1.0, Indices{index::a, index::e}, &F_ae_tilde, -1.0,
               Indices{index::n, index::m, index::f, index::a}, T_mn,
               Indices{index::n, index::m, index::f, index::e}, L_menf_t);

        // Tilde F_mi: canonical Eq. (40), local SI Eq. (S4); dimensions 4 * (m).
        std::array<Tensor<double, 1>, 4> F_mi_tilde_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            F_mi_tilde_list[i_idx] = F_mi_list[i_idx];

            Tensor<double, 3> T_mi_t("T_mi_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::n, index::e, index::f}, &T_mi_t, Indices{index::n, index::f, index::e}, T_mi_list[i_idx]);
            einsum(1.0, Indices{index::m}, &F_mi_tilde_list[i_idx], 1.0, Indices{index::m, index::n, index::e, index::f}, L_menf_t,
                    Indices{index::n, index::e, index::f}, T_mi_t);
        }
        L_menf_t = Tensor<double, 4>("released", 0, 0, 0, 0);

        // E_ei^{ma}: canonical Eq. (41), local SI Eq. (S5); dimensions 4 * (m, e, a).
        std::array<Tensor<double, 3>, 4> E_eima_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            // Bare-integral contribution in SI Eq. (S5): 2(me|ai) - (mi|ae).
            E_eima_list[i_idx] = g_meai_list[i_idx];
            E_eima_list[i_idx] *= 2.0;
            E_eima_list[i_idx] -= g_miae_list[i_idx];

            // Amplitude contribution in SI Eq. (S5): [2(me|nf) - (mf|ne)] U_ni^{fa}.
            einsum(1.0, Indices{index::m, index::e, index::a}, &E_eima_list[i_idx], 1.0, Indices{index::m, index::e, index::n, index::f}, L_menf,
                    Indices{index::n, index::f, index::a}, U_mi_list[i_idx]);
        }

        // F_ie^{ma}: canonical Eq. (42), local SI Eq. (S6); dimensions 4 * (m, e, a).
        std::array<Tensor<double, 3>, 4> F_iema_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            // Bare-integral contribution in SI Eq. (S6): (mi|ae).
            F_iema_list[i_idx] = g_miae_list[i_idx];

            // Amplitude contribution in SI Eq. (S6): -(mf|ne) T_ni^{af}.
            // Recover (mf|ne) from L_menf = 2.0 * (me|nf) - (mf|ne).
            Tensor<double, 3> T_mi_t("T_mi_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::n, index::f, index::a}, &T_mi_t, Indices{index::n, index::a, index::f}, T_mi_list[i_idx]);

            einsum(1.0, Indices{index::m, index::e, index::a}, &F_iema_list[i_idx], -2.0, Indices{index::m, index::e, index::n, index::f}, g_menf,
                    Indices{index::n, index::f, index::a}, T_mi_t);
            einsum(1.0, Indices{index::m, index::e, index::a}, &F_iema_list[i_idx], 1.0, Indices{index::m, index::e, index::n, index::f}, L_menf,
                    Indices{index::n, index::f, index::a}, T_mi_t);
        }

        // G_ij^{mn}: canonical Eq. (43), local SI Eq. (S7); dimensions 16 * (m, n).
        std::array<Tensor<double, 2>, 16> G_ijmn_list;
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                G_ijmn_list[i_idx * FOUR + j_idx] = Tensor<double, 2>("G_ijmn", nlmo_ijkl, nlmo_ijkl);
                // Bare-integral contribution in SI Eq. (S7): (mi|nj).
                einsum(0.0, Indices{index::m, index::n}, &G_ijmn_list[i_idx * FOUR + j_idx], 1.0, Indices{index::Q, index::m},
                        q_io_t1_list[i_idx], Indices{index::Q, index::n}, q_io_t1_list[j_idx]);

                // Amplitude contribution in SI Eq. (S7): (me|nf) T_ij^{ef}.
                einsum(1.0, Indices{index::m, index::n}, &G_ijmn_list[i_idx * FOUR + j_idx], 1.0, 
                        Indices{index::m, index::n, index::e, index::f}, g_menf_t, Indices{index::e, index::f}, T_ij_list[i_idx * FOUR + j_idx]);
            } // end j_idx
        } // end i_idx

        // H_ef^{ab}: canonical Eq. (44), local SI Eq. (S8); dimensions (e, f, a, b).
        Tensor<double, 4> H_efab("H_efab", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl); {
            // Bare-integral contribution in SI Eq. (S8): (ae|bf).
            Tensor<double, 4> H_efab_temp("H_efab_temp", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            einsum(0.0, Indices{index::e, index::a, index::f, index::b}, &H_efab_temp, 1.0, Indices{index::Q, index::e, index::a}, q_vv_t1, 
                    Indices{index::Q, index::f, index::b}, q_vv_t1);
            permute(Indices{index::e, index::f, index::a, index::b}, &H_efab, Indices{index::e, index::a, index::f, index::b}, H_efab_temp);

            // Amplitude contribution in SI Eq. (S8): (me|nf) T_mn^{ab}.
            einsum(1.0, Indices{index::e, index::f, index::a, index::b}, &H_efab, 1.0, 
                    Indices{index::m, index::n, index::e, index::f}, g_menf_t, Indices{index::m, index::n, index::a, index::b}, T_mn);
        }

        // I_eij^{mab}: canonical Eq. (45), local SI Eq. (S9).  Only the six
        // upper positional pairs enter Eq. (50).  Build the two raw orientations
        // for one pair, symmetrize immediately, and release the work tensors;
        // the former schedule retained sixteen raw rank-four tensors alongside
        // the final bank.
        std::array<Tensor<double, 4>, 16> I_eijmab_list;
        auto build_I_eijmab_raw = [&](int i_idx, int j_idx) {
            Tensor<double, 4> I_eijmab("I_eijmab raw", nlmo_ijkl, nqno_ijkl,
                                        nqno_ijkl, nqno_ijkl);

            // First term of SI Eq. (S9): [2(me|af) - (mf|ae)] T_ij^{fb}.
            einsum(0.0, Indices{index::m, index::e, index::a, index::b},
                   &I_eijmab, 1.0,
                   Indices{index::m, index::e, index::a, index::f}, L_mfae_t,
                   Indices{index::f, index::b},
                   T_ij_list[i_idx * FOUR + j_idx]);

            // Second term of SI Eq. (S9): -[2(me|ni) - (mi|ne)] T_nj^{ab}.
            einsum(1.0, Indices{index::m, index::e, index::a, index::b},
                   &I_eijmab, -1.0,
                   Indices{index::m, index::e, index::n}, L_menj_t_list[i_idx],
                   Indices{index::n, index::a, index::b}, T_mi_list[j_idx]);
            return I_eijmab;
        };
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            for (int j_idx = i_idx + 1; j_idx < FOUR; ++j_idx) {
                Tensor<double, 4> I_eijmab = build_I_eijmab_raw(i_idx, j_idx);
                Tensor<double, 4> I_ejimab = build_I_eijmab_raw(j_idx, i_idx);
                Tensor<double, 4> I_ejimba_buffer(
                    "I_ejimba", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                permute(Indices{index::m, index::e, index::b, index::a},
                        &I_ejimba_buffer,
                        Indices{index::m, index::e, index::a, index::b},
                        I_ejimab);

                const int ij_position = i_idx * FOUR + j_idx;
                I_eijmab_list[ij_position] = std::move(I_eijmab);
                I_eijmab_list[ij_position] += I_ejimba_buffer;

                // Final term of SI Eq. (S9): 1/4 [2(me|nf) - (mf|ne)]
                // Z_nij^{fab}; the factor includes P_ij^{ab}.
                Tensor<double, 4> Z_nij = form_Z_mij(i_idx, j_idx);
                einsum(1.0, Indices{index::m, index::e, index::a, index::b},
                       &I_eijmab_list[ij_position], 0.5,
                       Indices{index::m, index::e, index::n, index::f}, L_menf,
                       Indices{index::n, index::f, index::a, index::b}, Z_nij);
            }
        }

        // J_iej^{mab}: canonical Eq. (46), local SI Eq. (S10).  Build one
        // ordered-pair block and immediately consume both complementary (k,l)
        // orientations into a single rank-four residual.  This prevents the J
        // bank from coexisting with the I/K/L/M intermediates.
        Tensor<double, 4> R_J_ijkl("R_J_ijkl", nqno_ijkl, nqno_ijkl,
                                    nqno_ijkl, nqno_ijkl);
        R_J_ijkl.zero();
        {
            Tensor<double, 4> g_mfae_v("g_mfae_v", nlmo_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::e, index::a, index::f}, &g_mfae_v, Indices{index::m, index::f, index::e, index::a}, g_mfae);
            Tensor<double, 4> g_menf_v("g_menf_v", nlmo_ijkl, nqno_ijkl, nlmo_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::e, index::n, index::f}, &g_menf_v, Indices{index::m, index::f, index::n, index::e}, g_menf);

            auto build_J_iejmab = [&](int i_idx, int j_idx) {
                Tensor<double, 4> J_iejmab("J_iejmab", nlmo_ijkl, nqno_ijkl,
                                            nqno_ijkl, nqno_ijkl);

                // First term of SI Eq. (S10): (mf|ae) T_ij^{fb}.
                einsum(0.0, Indices{index::m, index::e, index::a, index::b},
                       &J_iejmab, 1.0,
                       Indices{index::m, index::e, index::a, index::f}, g_mfae_v,
                       Indices{index::f, index::b},
                       T_ij_list[i_idx * FOUR + j_idx]);

                // Second term of SI Eq. (S10): -(mi|ne) T_nj^{ab}.
                Tensor<double, 3> g_menj_u("g_menj_u", nlmo_ijkl, nqno_ijkl,
                                            nlmo_ijkl);
                permute(Indices{index::m, index::e, index::n}, &g_menj_u,
                        Indices{index::e, index::n, index::m},
                        g_menj_list[i_idx]);
                einsum(1.0, Indices{index::m, index::e, index::a, index::b},
                       &J_iejmab, -1.0,
                       Indices{index::m, index::e, index::n}, g_menj_u,
                       Indices{index::n, index::a, index::b}, T_mi_list[j_idx]);

                Tensor<double, 4> T_mij_t("T_mij_t", nlmo_ijkl, nqno_ijkl,
                                           nqno_ijkl, nqno_ijkl);
                auto T_mij = T_mij_view(i_idx, j_idx);
                permute(Indices{index::n, index::f, index::a, index::b},
                        &T_mij_t,
                        Indices{index::n, index::a, index::f, index::b}, T_mij);

                // Final term of SI Eq. (S10): -1/2 (mf|ne) T_nij^{afb}.
                einsum(1.0, Indices{index::m, index::e, index::a, index::b},
                       &J_iejmab, -0.5,
                       Indices{index::m, index::e, index::n, index::f}, g_menf_v,
                       Indices{index::n, index::f, index::a, index::b}, T_mij_t);
                return J_iejmab;
            };

            // Preserve the historical first representative for repeated
            // occupied indices, then schedule those representatives by their
            // ordered (i,j) source block.
            std::array<bool, 24> first_J_representative{};
            std::unordered_set<size_t> J_occupied_orbits;
            einsums::for_sequence<24UL>([&](auto perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] =
                    quadruple_permutations_[perm_idx];
                const size_t occupied_key = quadruplet_key(
                    ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                    ijkl_list[l_idx], naocc);
                first_J_representative[perm_idx] =
                    J_occupied_orbits.insert(occupied_key).second;
            });

            for (int source_i_idx = 0; source_i_idx < FOUR; ++source_i_idx) {
                for (int source_j_idx = 0; source_j_idx < FOUR; ++source_j_idx) {
                    if (source_i_idx == source_j_idx) continue;
                    Tensor<double, 4> J_iejmab =
                        build_J_iejmab(source_i_idx, source_j_idx);

                    einsums::for_sequence<24UL>([&](auto representative_perm_idx) {
                        auto &[i_idx, j_idx, k_idx, l_idx] =
                            quadruple_permutations_[representative_perm_idx];
                        if (!first_J_representative[representative_perm_idx] ||
                            i_idx != source_i_idx || j_idx != source_j_idx)
                            return;

                        const size_t occupied_key = quadruplet_key(
                            ijkl_list[i_idx], ijkl_list[j_idx], ijkl_list[k_idx],
                            ijkl_list[l_idx], naocc);
                        Tensor<double, 4> R_ijkl_perm(
                            "R_ijkl_perm J", nqno_ijkl, nqno_ijkl,
                            nqno_ijkl, nqno_ijkl);
                        R_ijkl_perm.zero();
                        Tensor<double, 4> R_ijkl_buffer_a(
                            "R_ijkl_buffer_a J", nqno_ijkl, nqno_ijkl,
                            nqno_ijkl, nqno_ijkl);
                        Tensor<double, 4> R_ijkl_buffer_b(
                            "R_ijkl_buffer_b J", nqno_ijkl, nqno_ijkl,
                            nqno_ijkl, nqno_ijkl);

                        Tensor<double, 4> T_mkl_t(
                            "T_mkl_t", nlmo_ijkl, nqno_ijkl, nqno_ijkl,
                            nqno_ijkl);
                        auto T_mkl = T_mij_view(k_idx, l_idx);
                        permute(Indices{index::m, index::e, index::c, index::d},
                                &T_mkl_t,
                                Indices{index::m, index::c, index::e, index::d},
                                T_mkl);

                        // J_iej^{mab} contribution to canonical Eq. (50), local
                        // SI Eq. (S31).
                        einsum(0.0,
                               Indices{index::a, index::b, index::c, index::d},
                               &R_ijkl_buffer_a, 1.0,
                               Indices{index::m, index::e, index::a, index::b},
                               J_iejmab,
                               Indices{index::m, index::e, index::c, index::d},
                               T_mkl_t);
                        permute(Indices{index::c, index::b, index::a, index::d},
                                &R_ijkl_buffer_b,
                                Indices{index::a, index::b, index::c, index::d},
                                R_ijkl_buffer_a);
                        R_ijkl_buffer_a *= -0.5;
                        R_ijkl_perm += R_ijkl_buffer_a;
                        R_ijkl_buffer_b *= -1.0;
                        R_ijkl_perm += R_ijkl_buffer_b;

                        einsums::for_sequence<24UL>([&](auto target_perm_idx) {
                            auto &[target_i_idx, target_j_idx, target_k_idx,
                                   target_l_idx] =
                                quadruple_permutations_[target_perm_idx];
                            const size_t target_key = quadruplet_key(
                                ijkl_list[target_i_idx], ijkl_list[target_j_idx],
                                ijkl_list[target_k_idx], ijkl_list[target_l_idx],
                                naocc);
                            if (target_key != occupied_key) return;
                            permute(
                                Indices{index::a, index::b, index::c, index::d},
                                &R_ijkl_buffer_a,
                                std::get<target_perm_idx>(einsum_indices),
                                R_ijkl_perm);
                            R_J_ijkl += R_ijkl_buffer_a;
                        });
                    });
                }
            }
        }

        std::vector<int> needed_triplet_permutations;
        needed_triplet_permutations.reserve(ijk_to_ijkl_perm_idx.size());
        for (const auto& entry : ijk_to_ijkl_perm_idx) {
            needed_triplet_permutations.push_back(
                static_cast<int>(entry.second));
        }
        std::sort(needed_triplet_permutations.begin(),
                  needed_triplet_permutations.end());
        needed_triplet_permutations.erase(
            std::unique(needed_triplet_permutations.begin(),
                        needed_triplet_permutations.end()),
            needed_triplet_permutations.end());

        // K_ijk^{amn}: canonical Eq. (47), local SI Eq. (S11).  Only final
        // representatives addressable through ijk_to_ijkl_perm_idx are retained.
        // Raw tensors are cached until their final symmetrization use and then
        // released, instead of keeping raw and final 24-banks simultaneously.
        std::unordered_map<int, Tensor<double, 3>> K_ijkamn_map;
        {
            auto K_twin = [&](int perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] =
                    quadruple_permutations_[perm_idx];
                return static_cast<int>(ijk_to_ijkl_perm_idx.at(triplet_key(
                    ijkl_list[i_idx], ijkl_list[k_idx], ijkl_list[j_idx],
                    naocc)));
            };
            auto build_K_raw = [&](int perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] =
                    quadruple_permutations_[perm_idx];
                Tensor<double, 3> K_ijkamn("K_ijkamn raw", nqno_ijkl,
                                            nlmo_ijkl, nlmo_ijkl);

                // First term of SI Eq. (S11): (me|nk) T_ij^{ae}.
                einsum(0.0, Indices{index::a, index::m, index::n},
                       &K_ijkamn, 1.0, Indices{index::a, index::e},
                       T_ij_list[i_idx * FOUR + j_idx],
                       Indices{index::e, index::m, index::n},
                       g_menj_list[k_idx]);

                // Second term of SI Eq. (S11): 1/2 T_ijk^{aef} (me|nf).
                // T_mni_list[k](i,j,...) is the already projected T_ijk block.
                auto T_ijk = T_mni_list[k_idx](
                    occupied_domain_positions[i_idx],
                    occupied_domain_positions[j_idx], All, All, All);
                Tensor<double, 3> K_second_mna("K second (m,n,a)", nlmo_ijkl,
                                                nlmo_ijkl, nqno_ijkl);
                einsum(0.0, Indices{index::m, index::n, index::a},
                       &K_second_mna, 0.5,
                       Indices{index::m, index::n, index::e, index::f},
                       g_menf_t, Indices{index::a, index::e, index::f}, T_ijk);
                Tensor<double, 3> K_second_amn("K second (a,m,n)", nqno_ijkl,
                                                nlmo_ijkl, nlmo_ijkl);
                permute(Indices{index::a, index::m, index::n}, &K_second_amn,
                        Indices{index::m, index::n, index::a}, K_second_mna);
                K_ijkamn += K_second_amn;
                return K_ijkamn;
            };

            std::unordered_map<int, int> raw_use_count;
            for (int perm_idx : needed_triplet_permutations) {
                ++raw_use_count[perm_idx];
                ++raw_use_count[K_twin(perm_idx)];
            }
            std::unordered_map<int, Tensor<double, 3>> raw_cache;
            auto get_K_raw = [&](int perm_idx) -> Tensor<double, 3>& {
                auto found = raw_cache.find(perm_idx);
                if (found == raw_cache.end()) {
                    found = raw_cache.emplace(perm_idx, build_K_raw(perm_idx)).first;
                }
                return found->second;
            };
            auto release_K_raw = [&](int perm_idx) {
                if (--raw_use_count[perm_idx] == 0) raw_cache.erase(perm_idx);
            };

            for (int perm_idx : needed_triplet_permutations) {
                const int twin_perm_idx = K_twin(perm_idx);
                Tensor<double, 3> K_ijkamn = get_K_raw(perm_idx);
                Tensor<double, 3> K_ikjanm("K_ikjanm", nqno_ijkl, nlmo_ijkl,
                                            nlmo_ijkl);
                permute(Indices{index::a, index::n, index::m}, &K_ikjanm,
                        Indices{index::a, index::m, index::n},
                        get_K_raw(twin_perm_idx));
                K_ijkamn += K_ikjanm;
                K_ijkamn_map.emplace(perm_idx, std::move(K_ijkamn));
                release_K_raw(perm_idx);
                release_K_raw(twin_perm_idx);
            }
        }

        // L_ijk^{abm}: canonical Eq. (48), local SI Eqs. (S12)--(S13).
        // Apply the same reference-counted orbit schedule as K.
        std::unordered_map<int, Tensor<double, 3>> L_ijkabm_map;
        {
            auto L_twin = [&](int perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] =
                    quadruple_permutations_[perm_idx];
                return static_cast<int>(ijk_to_ijkl_perm_idx.at(triplet_key(
                    ijkl_list[j_idx], ijkl_list[i_idx], ijkl_list[k_idx],
                    naocc)));
            };
            auto build_L_raw = [&](int perm_idx) {
                auto &[i_idx, j_idx, k_idx, l_idx] =
                    quadruple_permutations_[perm_idx];
                const int i = ijkl_list[i_idx];
                const int j = ijkl_list[j_idx];
                const int k = ijkl_list[k_idx];
                const int ij = i_j_to_ij_[i][j];

                Tensor<double, 3> F_mae("F_mae", nlmo_ijkl, nqno_ijkl,
                                         nqno_ijkl);
                permute(Indices{index::m, index::a, index::e}, &F_mae,
                        Indices{index::m, index::e, index::a},
                        F_iema_list[k_idx]);
                Tensor<double, 3> EF_mea_sum = E_eima_list[i_idx];
                EF_mea_sum += F_iema_list[i_idx];
                Tensor<double, 3> EF_mae("EF_mae", nlmo_ijkl, nqno_ijkl,
                                          nqno_ijkl);
                permute(Indices{index::m, index::a, index::e}, &EF_mae,
                        Indices{index::m, index::e, index::a}, EF_mea_sum);

                Tensor<double, 3> L_ijkabm("L_ijkabm raw", nlmo_ijkl,
                                            nqno_ijkl, nqno_ijkl);
                L_ijkabm.zero();

                // First term of SI Eq. (S12): (mf|ae) T_ijk^{ebf}.
                // Reuse T_mni_list[k](i,j,...) instead of projecting T3 again.
                auto T_ijk = T_mni_list[k_idx](
                    occupied_domain_positions[i_idx],
                    occupied_domain_positions[j_idx], All, All, All);
                Tensor<double, 3> T_ijk_t("T_ijk_t", nqno_ijkl, nqno_ijkl,
                                           nqno_ijkl);
                permute(Indices{index::f, index::e, index::b}, &T_ijk_t,
                        Indices{index::e, index::b, index::f}, T_ijk);
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm,
                       1.0, Indices{index::m, index::a, index::f, index::e},
                       g_mfae_u, Indices{index::f, index::e, index::b}, T_ijk_t);

                // Second term of SI Eq. (S12):
                // 1/2 (E_ei^{ma} + F_ie^{ma}) T_jk^{be}.
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm,
                       0.5, Indices{index::m, index::a, index::e}, EF_mae,
                       Indices{index::b, index::e},
                       T_ij_list[j_idx * FOUR + k_idx]);

                // Third term of SI Eq. (S12): F_ke^{ma} T_ij^{eb}.
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm,
                       1.0, Indices{index::m, index::a, index::e}, F_mae,
                       Indices{index::e, index::b},
                       T_ij_list[i_idx * FOUR + j_idx]);

                // Fourth term of SI Eq. (S12): -1/2 G_ki^{mn} T_nj^{ab}.
                einsum(1.0, Indices{index::m, index::a, index::b}, &L_ijkabm,
                       -0.5, Indices{index::m, index::n},
                       G_ijmn_list[k_idx * FOUR + i_idx],
                       Indices{index::n, index::a, index::b}, T_mi_list[j_idx]);

                // Pair-domain delta L term, SI Eq. (S13) and main-text
                // Eqs. (98)--(99).  Each two-index projection is a GEMM.
                auto S_ijkl_ij = submatrix_rows_and_cols(
                    *S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[ij]);
                S_ijkl_ij = linalg::triplet(X_qno_[ijkl], S_ijkl_ij,
                                            X_pno_[ij], true, false, false);
                Tensor<double, 2> S_ijkl_ij_ein("S_ijkl_ij_ein", nqno_ijkl,
                                                 n_pno_[ij]);
                ::memcpy(S_ijkl_ij_ein.data(), S_ijkl_ij->get_pointer(),
                         static_cast<size_t>(nqno_ijkl) * n_pno_[ij] *
                             sizeof(double));

                const int k_ij = lmopair_to_lmos_dense_[ij][k];
                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    const int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    const int m_ij = lmopair_to_lmos_dense_[ij][m];
                    if (m_ij == -1) continue;

                    auto delta_L_pair_slice =
                        delta_L_ijk_abm_list[ij](k_ij, m_ij, All, All);
                    Tensor<double, 2> delta_L_qno_pno(
                        "delta_L_qno_pno", nqno_ijkl, n_pno_[ij]);
                    einsum(0.0, Indices{index::c, index::b},
                           &delta_L_qno_pno, 1.0,
                           Indices{index::c, index::a}, S_ijkl_ij_ein,
                           Indices{index::a, index::b}, delta_L_pair_slice);
                    auto L_ijkabm_slice = L_ijkabm(m_ijkl, All, All);
                    einsum(1.0, Indices{index::a, index::d},
                           &L_ijkabm_slice, 1.0,
                           Indices{index::a, index::b}, delta_L_qno_pno,
                           Indices{index::d, index::b}, S_ijkl_ij_ein);
                }
                return L_ijkabm;
            };

            std::unordered_map<int, int> raw_use_count;
            for (int perm_idx : needed_triplet_permutations) {
                ++raw_use_count[perm_idx];
                ++raw_use_count[L_twin(perm_idx)];
            }
            std::unordered_map<int, Tensor<double, 3>> raw_cache;
            auto get_L_raw = [&](int perm_idx) -> Tensor<double, 3>& {
                auto found = raw_cache.find(perm_idx);
                if (found == raw_cache.end()) {
                    found = raw_cache.emplace(perm_idx, build_L_raw(perm_idx)).first;
                }
                return found->second;
            };
            auto release_L_raw = [&](int perm_idx) {
                if (--raw_use_count[perm_idx] == 0) raw_cache.erase(perm_idx);
            };

            for (int perm_idx : needed_triplet_permutations) {
                const int twin_perm_idx = L_twin(perm_idx);
                Tensor<double, 3> L_ijkabm = get_L_raw(perm_idx);
                Tensor<double, 3> L_jikbam("L_jikbam", nlmo_ijkl, nqno_ijkl,
                                            nqno_ijkl);
                permute(Indices{index::m, index::a, index::b}, &L_jikbam,
                        Indices{index::m, index::b, index::a},
                        get_L_raw(twin_perm_idx));
                L_ijkabm += L_jikbam;
                L_ijkabm_map.emplace(perm_idx, std::move(L_ijkabm));
                release_L_raw(perm_idx);
                release_L_raw(twin_perm_idx);
            }
        }

        // M_ejk^{abc}: canonical Eq. (49), local SI Eqs. (S14)--(S15).
        // Eq. (50) consumes only the six upper positional pairs.  Construct the
        // two raw orientations for one pair and symmetrize directly, avoiding
        // the former 16-tensor raw bank and ten unused final entries.
        std::array<Tensor<double, 4>, 16> M_ejkabc_list;
        {
            Tensor<double, 4> H_efab_transpose("H_efab_transpose", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::e, index::a, index::b, index::f}, &H_efab_transpose, Indices{index::e, index::f, index::a, index::b}, H_efab);

            auto build_M_raw = [&](int j_idx, int k_idx) {
                const int j = ijkl_list[j_idx];
                const int k = ijkl_list[k_idx];
                const int jk = i_j_to_ij_[j][k];
                Tensor<double, 4> M_ejkabc("M_ejkabc raw", nqno_ijkl,
                                            nqno_ijkl, nqno_ijkl, nqno_ijkl);

                // First term of SI Eq. (S14): 1/2 H_ef^{ab} T_jk^{fc}.
                einsum(0.0, Indices{index::e, index::a, index::b, index::c},
                       &M_ejkabc, 0.5,
                       Indices{index::e, index::a, index::b, index::f},
                       H_efab_transpose, Indices{index::f, index::c},
                       T_ij_list[j_idx * FOUR + k_idx]);

                // Pair-domain delta M term, SI Eq. (S15) and main-text
                // Eqs. (101)--(102).
                auto S_ijkl_jk = submatrix_rows_and_cols(
                    *S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[jk]);
                S_ijkl_jk = linalg::triplet(X_qno_[ijkl], S_ijkl_jk,
                                            X_pno_[jk], true, false, false);
                M_ejkabc += matmul_4d(delta_M_ejk_abc_list[jk], S_ijkl_jk,
                                      n_pno_[jk], n_qno_[ijkl]);
                return M_ejkabc;
            };

            for (int j_idx = 0; j_idx < FOUR; ++j_idx) {
                for (int k_idx = j_idx + 1; k_idx < FOUR; ++k_idx) {
                    Tensor<double, 4> M_ejkabc = build_M_raw(j_idx, k_idx);
                    Tensor<double, 4> M_ekjacb_source = build_M_raw(k_idx, j_idx);
                    Tensor<double, 4> M_ekjacb("M_ekjacb", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
                    permute(Indices{index::e, index::a, index::b, index::c},
                            &M_ekjacb,
                            Indices{index::e, index::a, index::c, index::b},
                            M_ekjacb_source);
                    const int jk_position = j_idx * FOUR + k_idx;
                    M_ejkabc_list[jk_position] = std::move(M_ejkabc);
                    M_ejkabc_list[jk_position] += M_ekjacb;
                }
            }
        }

        // => Form all possible R_ijkl's over unique (i, j, k, l)
        Tensor<double, 4> R_ijkl("R_ijkl", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
        R_ijkl.zero();

        // Terms with (i, jkl)-type symmetry
        for (int i_idx = 0; i_idx < FOUR; ++i_idx) {
            int i = ijkl_list[i_idx];
            auto &[j, k, l] = complementary_triplets[i_idx];

            Tensor<double, 4> R_ijkl_buffer("R_ijkl_buffer", nqno_ijkl, nqno_ijkl, nqno_ijkl, nqno_ijkl);
            R_ijkl_buffer.zero();

            // Tilde-F_ae contribution to canonical Eq. (50), arranged as local
            // SI Eqs. (S22)--(S23).
            // (permutationally adapted coefficient is +1)
            Tensor<double, 4> T_ijkl = quadruples_permuter(T_iajbkcld_[ijkl], i, j, k, l);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer, 1.0,
                    Indices{index::a, index::e}, F_ae_tilde, Indices{index::e, index::b, index::c, index::d}, T_ijkl);

            // Tilde-F_mi contribution to canonical Eq. (50), arranged as local
            // SI Eqs. (S24)--(S25); spin/permutation adaptation gives coefficient -1.
            // The raw, occupied-permuted, and alpha complement banks formerly
            // occupied 12*nlmo*nqno^4 doubles.  Retain only the raw four-bank
            // source and consume one sorted/alpha m slice at a time.
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                const int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                Tensor<double, 4> T_mijk_unsorted =
                    T_nijk_complement_unsorted[i_idx](m_ijkl, All, All, All, All);
                Tensor<double, 4> T_mijk =
                    quadruples_permuter(T_mijk_unsorted, m, j, k, l);

                const size_t nqno2 = static_cast<size_t>(nqno_ijkl) * nqno_ijkl;
                C_DAXPY(nqno2 * nqno2, -F_mi_tilde_list[i_idx](m_ijkl),
                        T_mijk.data(), 1, R_ijkl_buffer.data(), 1);

                Tensor<double, 4> alpha_mijk = form_alpha_ijkl(T_mijk);
                auto E_eia = E_eima_list[i_idx](m_ijkl, All, All);
                einsum(1.0, Indices{index::a, index::b, index::c, index::d},
                       &R_ijkl_buffer, 0.5,
                       Indices{index::e, index::a}, E_eia,
                       Indices{index::e, index::b, index::c, index::d},
                       alpha_mijk);
            }

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
            Tensor<double, 4> R_ijkl_buffer_c("R_ijkl_buffer_c", n_xpno_[kl], n_xpno_[kl], n_xpno_[kl], n_xpno_[kl]);

            R_ijkl_buffer_a.zero();
            R_ijkl_buffer_b.zero();
            R_ijkl_buffer_c.zero();

            // A_ej^{ab} contribution to canonical Eq. (50), local SI Eq. (S20).
            const size_t ikl_dense = triplet_key(i, k, l, naocc);
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

            const size_t jkl_dense = triplet_key(j, k, l, naocc);
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

            // B_ij^{am} contribution to canonical Eq. (50), local SI Eq. (S21).
            // Fixing the second occupied axis of T_mni leaves a strided view in
            // the leading m dimension. Pack this one block so both B*T3 calls
            // are ordinary GEMMs over a contiguous (m,b,c,d) operand.
            auto T_mkl_view = T_mij_view(k_idx, l_idx);
            Tensor<double, 4> T_mkl("T_mkl contiguous", nlmo_ijkl,
                                     nqno_ijkl, nqno_ijkl, nqno_ijkl);
            permute(Indices{index::m, index::b, index::c, index::d}, &T_mkl,
                    Indices{index::m, index::b, index::c, index::d},
                    T_mkl_view);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, -1.0,
                    Indices{index::m, index::a}, B_ijam_list[i_idx * FOUR + j_idx], Indices{index::m, index::b, index::c, index::d}, T_mkl);

            einsum(0.0, Indices{index::b, index::a, index::c, index::d}, &R_ijkl_buffer_b, -1.0,
                    Indices{index::m, index::b}, B_ijam_list[j_idx * FOUR + i_idx], Indices{index::m, index::a, index::c, index::d}, T_mkl);

            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 1, 0, 2, 3);

            // F_ie^{ma} contribution to canonical Eq. (50), local SI Eq. (S27).
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                // Order both (m, i, j, l) and (m, j, k, l) for the permutation lookup.
                Tensor<double, 4> T_mikl = T_nijk_complement_unsorted[j_idx](m_ijkl, All, All, All, All);
                Tensor<double, 4> T_mjkl = T_nijk_complement_unsorted[i_idx](m_ijkl, All, All, All, All);

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
            
            // G_ij^{mn} contribution to canonical Eq. (50), local SI Eq. (S28).
            // The XPNO projection follows main-text Eqs. (83)--(85).
            for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                
                for (int n_ijkl = 0; n_ijkl < nlmo_ijkl; ++n_ijkl) {
                    int n = lmoquadruplet_to_lmos_[ijkl][n_ijkl];
                    int mn = i_j_to_ij_[m][n];
                    const size_t mnkl_idx = quadruplet_key(m, n, k, l, naocc);
                    int mnkl = i_j_k_l_to_ijkl_.count(mnkl_idx) ? i_j_k_l_to_ijkl_[mnkl_idx] : -1;
                    if (mn == -1 || mnkl == -1) continue;
                    
                    Tensor<double, 4> T_mnkl;
                    if (stream_xpno_t4_) {
                        // Low-memory path: project only the block being consumed.
                        // Perform four canonical GEMMs and apply both occupied
                        // pair orientations only in the smaller XPNO target
                        // space, so no persistent [kl][mn] bank or old-space
                        // rank-four permutation is required.
                        const int canonical_kl = (k <= l) ? kl : ij_to_ji_[kl];
                        auto S_mnkl_kl = submatrix_rows_and_cols(
                            *S_pao_, lmoquadruplet_to_paos_[mnkl],
                            lmopair_to_paos_ext_[canonical_kl]);
                        S_mnkl_kl = linalg::triplet(
                            X_qno_[mnkl], S_mnkl_kl, X_xpno_[canonical_kl],
                            true, false, false);
                        T_mnkl = matmul_4d_permuted(
                            T_iajbkcld_[mnkl], S_mnkl_kl->transpose(),
                            n_qno_[mnkl], n_xpno_[canonical_kl], m, n, k, l);
                    } else {
                        T_mnkl = Tensor<double, 4>(
                            "T_mnkl", n_xpno_[kl], n_xpno_[kl], n_xpno_[kl],
                            n_xpno_[kl]);
                        int lk = ij_to_ji_[kl], nm = ij_to_ji_[mn];

                        if (m > n && k > l) {
                            permute(Indices{index::b, index::a, index::d, index::c},
                                    &T_mnkl,
                                    Indices{index::a, index::b, index::c, index::d},
                                    T_mnkl_xpno_[lk][nm]);
                        } else if (m > n) {
                            permute(Indices{index::b, index::a, index::c, index::d},
                                    &T_mnkl,
                                    Indices{index::a, index::b, index::c, index::d},
                                    T_mnkl_xpno_[kl][nm]);
                        } else if (k > l) {
                            permute(Indices{index::a, index::b, index::d, index::c},
                                    &T_mnkl,
                                    Indices{index::a, index::b, index::c, index::d},
                                    T_mnkl_xpno_[lk][mn]);
                        } else {
                            T_mnkl = T_mnkl_xpno_[kl][mn];
                        }
                    }

                    const size_t nxpno = static_cast<size_t>(n_xpno_[kl]);
                    const size_t nxpno2 = nxpno * nxpno;
                    const size_t length = nxpno2 * nxpno2;
                    C_DAXPY(length, (G_ijmn_list[i_idx * FOUR + j_idx])(m_ijkl, n_ijkl), T_mnkl.data(), 1, R_ijkl_buffer_c.data(), 1);
                } // end n_ijkl
            } // end m_ijkl
            // Flush G contributions
            auto S_ijkl_kl = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_ext_[kl]);
            S_ijkl_kl = linalg::triplet(X_qno_[ijkl], S_ijkl_kl, X_xpno_[kl], true, false, false);
            R_ijkl_buffer_a += matmul_4d(R_ijkl_buffer_c, S_ijkl_kl, n_xpno_[kl], n_qno_[ijkl]);

            // H_ef^{ab} contribution to canonical Eq. (50), local SI Eq. (S29).
            Tensor<double, 4> T_ijkl = quadruples_permuter(T_iajbkcld_[ijkl], i, j, k, l);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, Indices{index::e, index::f, index::a, index::b}, H_efab,
                    Indices{index::e, index::f, index::c, index::d}, T_ijkl);

            // I_eij^{mab} contribution to canonical Eq. (50), local SI Eq. (S30).
            Tensor<double, 4> Z_mkl = form_Z_mij(k_idx, l_idx);
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 0.5, Indices{index::m, index::e, index::a, index::b}, 
                    I_eijmab_list[i_idx * FOUR + j_idx], Indices{index::m, index::e, index::c, index::d}, Z_mkl);

            // K_ijk^{amn} contribution to canonical Eq. (50), local SI Eq. (S32).
            // (becomes) T_{mnk}^{abc} K_{lij}^{dmn}
            size_t lij_perm_idx =
                ijk_to_ijkl_perm_idx.at(triplet_key(l, i, j, naocc));
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, 1.0, Indices{index::m, index::n, index::a, index::b, index::c}, 
                    T_mni_list[k_idx], Indices{index::d, index::m, index::n},
                    K_ijkamn_map.at(static_cast<int>(lij_perm_idx)));

            size_t kij_perm_idx =
                ijk_to_ijkl_perm_idx.at(triplet_key(k, i, j, naocc));
            einsum(0.0, Indices{index::a, index::b, index::d, index::c}, &R_ijkl_buffer_b, 1.0, Indices{index::m, index::n, index::a, index::b, index::d},
                    T_mni_list[l_idx], Indices{index::c, index::m, index::n},
                    K_ijkamn_map.at(static_cast<int>(kij_perm_idx)));
            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 0, 1, 3, 2);

            // L_ijk^{abm} contribution to canonical Eq. (50), local SI Eq. (S33).
            size_t ijk_perm_idx =
                ijk_to_ijkl_perm_idx.at(triplet_key(i, j, k, naocc));
            einsum(1.0, Indices{index::a, index::b, index::c, index::d}, &R_ijkl_buffer_a, -1.0, Indices{index::m, index::a, index::b}, 
                    L_ijkabm_map.at(static_cast<int>(ijk_perm_idx)),
                    Indices{index::m, index::c, index::d}, T_mi_list[l_idx]);
            
            size_t ijl_perm_idx =
                ijk_to_ijkl_perm_idx.at(triplet_key(i, j, l, naocc));
            einsum(0.0, Indices{index::a, index::b, index::d, index::c}, &R_ijkl_buffer_b, -1.0, Indices{index::m, index::a, index::b}, 
                    L_ijkabm_map.at(static_cast<int>(ijl_perm_idx)),
                    Indices{index::m, index::d, index::c}, T_mi_list[k_idx]);
            R_ijkl_buffer_a += quadruples_permuter(R_ijkl_buffer_b, 0, 1, 3, 2);

            // M_ejk^{abc} contribution to canonical Eq. (50), local SI Eq. (S34).
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

        // The J contribution was produced pair-by-pair before the K/L/M banks;
        // add it in the historical final position to preserve accumulation order.
        R_ijkl += R_J_ijkl;

        if (disk_qno_integrals_) {
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
    T_m_ijkl_.resize(n_lmo_quadruplets);

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

            // Project T1 into the pair, triplet, and quadruplet virtual spaces.
    #pragma omp parallel for schedule(dynamic, 1)
            for (int ij = 0; ij < n_lmo_pairs; ++ij) {
                auto &[i, j] = ij_to_i_j_[ij];

                int nlmo_ij = lmopair_to_lmos_[ij].size();
                int npno_ij = n_pno_[ij];
                int ij_idx = (i <= j) ? ij : ij_to_ji_[ij];
                
                T_n_ij_[ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);

                for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                    int n = lmopair_to_lmos_[ij][n_ij];
                    auto T_n_temp = linalg::doublet(S_pno_ij_nn_[ij_idx][n], T_ia_[n], false, false);
                    
                    for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                        (*T_n_ij_[ij])(n_ij, a_ij) = (*T_n_temp)(a_ij, 0);
                    } // end a_ij
                } // end n_ij
            }

    #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];
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
                int nlmo_ijkl = lmoquadruplet_to_lmos_[ijkl].size();

                T_m_ijkl_[ijkl] = Tensor<double, 2>("T_m_ijkl", nlmo_ijkl, n_qno_[ijkl]);
                
                for (int m_ijkl = 0; m_ijkl < nlmo_ijkl; ++m_ijkl) {
                    int m = lmoquadruplet_to_lmos_[ijkl][m_ijkl];
                    int mm = i_j_to_ij_[m][m];

                    auto S_ijkl_mm = submatrix_rows_and_cols(*S_pao_, lmoquadruplet_to_paos_[ijkl], lmopair_to_paos_[mm]);
                    S_ijkl_mm = linalg::triplet(X_qno_[ijkl], S_ijkl_mm, X_pno_[mm], true, false, false);

                    auto T_m_temp = linalg::doublet(S_ijkl_mm, T_ia_[m]);

                    for (int a_ijkl = 0; a_ijkl < n_qno_[ijkl]; ++a_ijkl) {
                        (T_m_ijkl_[ijkl])(m_ijkl, a_ijkl) = (*T_m_temp)(a_ijkl, 0);
                    }
                } // end m_ijkl
            } // end ijkl
            
            // Create T_iajbkc_clone intermediate
    #pragma omp parallel for schedule(dynamic, 1)
            for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                int ijk = sorted_triplets_[ijk_sorted];

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

                T_iajbkc_clone_[ijk] = triples_spin_desummation(triples_spin_summation(T_iajbkc_clone_[ijk]));
                U_iajbkc_[ijk] = triples_spin_summation(T_iajbkc_clone_[ijk]);
                ::memcpy(T_iajbkc_[ijk]->get_pointer(), T_iajbkc_clone_[ijk].data(), n_tno_[ijk] * n_tno_[ijk] * n_tno_[ijk] * sizeof(double));
            }

            // spin-adapt and then de-adapt quadruples amplitudes
        #pragma omp parallel for schedule(dynamic, 1)
            for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                int ijkl = sorted_quadruplets_[ijkl_sorted];

                T_iajbkcld_[ijkl] = quadruples_spin_desummation(quadruples_spin_summation(T_iajbkcld_[ijkl]));
            }

            // compute quadruples amplitude
            timer_on("DLPNO-CCSDTQ : R_iajbkcld");
            if (miter == 0) {
                if (!stream_xpno_t4_) form_T_mnkl_xpno();
                compute_quadruples_residual(R_iajbkcld);
                // Equation (85) consumes the XPNO-projected amplitudes completely;
                // release them before the T3/T2/T1 residual phases begin.
                if (!stream_xpno_t4_) T_mnkl_xpno_.clear();

            #pragma omp parallel for schedule(dynamic, 1)
                for (int ijkl_sorted = 0; ijkl_sorted < n_lmo_quadruplets; ++ijkl_sorted) {
                    int ijkl = sorted_quadruplets_[ijkl_sorted];

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

                    const double residual_rms = rmsd(R_iajbkcld[ijkl], zero_tensor);
                    double alpha = (fabs(residual_rms) > fabs(R_iajbkcld_rms[ijkl]))
                                       ? quadruples_damping_ratio_
                                       : 0.0;

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

                    R_iajbkcld_rms[ijkl] = residual_rms;
                    r_curr4 += R_iajbkcld_rms[ijkl] * R_iajbkcld_rms[ijkl];
                }
                r_curr4 = n_lmo_quadruplets > 0
                              ? std::sqrt(r_curr4 / n_lmo_quadruplets)
                              : 0.0;
            }
            timer_off("DLPNO-CCSDTQ : R_iajbkcld");

            // compute triples amplitude
            timer_on("DLPNO-CCSDTQ : R_iajbkc");
            if (miter == 0) {
                // form_T_mnk();
                add_t4_to_triples_residual(R_iajbkc);

                // spin adapt and then de-adapt triples residual
        #pragma omp parallel for schedule(dynamic, 1)
                for (int ijk_sorted = 0; ijk_sorted < n_lmo_triplets; ++ijk_sorted) {
                    int ijk = sorted_triplets_[ijk_sorted];

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
                    double alpha = (fabs(R_iajbkc[ijk]->rms()) > fabs(R_iajbkc_rms[ijk])) ? quadruples_damping_ratio_ : 0.0;

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
            add_t4_to_doubles_residual(R_iajb, Rn_iajb, R_iajb_buffer);

            // Update doubles amplitude
            r_curr2 = 0.0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+ : r_curr2)
            for (int ij = 0; ij < n_lmo_pairs; ++ij) {
                auto &[i, j] = ij_to_i_j_[ij];
                double alpha = (fabs(R_iajb[ij]->rms()) > fabs(R_iajb_rms[ij])) ? quadruples_damping_ratio_ : 0.0;

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
                double alpha = (fabs(R_ia[i]->rms()) > fabs(R_ia_rms[i])) ? quadruples_damping_ratio_ : 0.0;

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
        const size_t nelements = T_ia_.size() + T_iajb_.size() + T_iajbkc_.size();

        std::vector<SharedMatrix> T_vecs;
        T_vecs.reserve(nelements);
        T_vecs.insert(T_vecs.end(), T_ia_.begin(), T_ia_.end());
        T_vecs.insert(T_vecs.end(), T_iajb_.begin(), T_iajb_.end());
        T_vecs.insert(T_vecs.end(), T_iajbkc_.begin(), T_iajbkc_.end());

        std::vector<SharedMatrix> R_vecs;
        R_vecs.reserve(nelements);
        R_vecs.insert(R_vecs.end(), R_ia.begin(), R_ia.end());
        R_vecs.insert(R_vecs.end(), R_iajb.begin(), R_iajb.end());
        R_vecs.insert(R_vecs.end(), R_iajbkc.begin(), R_iajbkc.end());

        // Flatten the native Einsums T4/R4 tensors directly into the two vectors
        // required by DIIS. This avoids persistent duplicate Psi4 Matrix copies.
        auto T_vecs_flat = flatten_ccsdtq_diis(T_vecs, T_iajbkcld_, extrapolate_t4_);
        auto R_vecs_flat = flatten_ccsdtq_diis(R_vecs, R_iajbkcld, extrapolate_t4_);

        if (iteration == 1) {
            diis.set_error_vector_size(R_vecs_flat);
            diis.set_vector_size(T_vecs_flat);
        }

        diis.add_entry(R_vecs_flat.get(), T_vecs_flat.get());
        diis.extrapolate(T_vecs_flat.get());

        copy_ccsdtq_diis(T_vecs_flat, T_vecs, T_iajbkcld_, extrapolate_t4_);

        // evaluate energy and convergence
        e_prev = e_curr;
        e_curr = 0.0;
        e_weak = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr, e_weak)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];
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
    DLPNOCCSDT_Q::compute_energy();

    timer_on("DLPNO-CCSDTQ");

    einsums::profile::initialize();

    print_header();

    // Build the XPNOs used to evaluate the projected T_mnkl contractions of
    // main-text Eqs. (83)--(85).
    double xpno_tolerance = options_.get_double("T_CUT_XPNO");
    xpno_transform(xpno_tolerance);

    timer_on("DLPNO-CCSDTQ : Estimate Memory");
    estimate_memory();
    timer_off("DLPNO-CCSDTQ : Estimate Memory");

    // estimate_memory() may activate disk-backed QNO integrals, so defer all
    // full-quadruples PSIO setup until the effective storage policy is known.
    if (write_qia_pno_) {
        psio_->open(PSIF_DLPNO_QIA_PNO, PSIO_OPEN_OLD);
    }
    if (write_qab_pno_) {
        psio_->open(PSIF_DLPNO_QAB_PNO, PSIO_OPEN_OLD);
    }
    if (disk_ints_) {
        psio_->open(PSIF_DLPNO_QIA_TNO, PSIO_OPEN_OLD);
        psio_->open(PSIF_DLPNO_QAB_TNO, PSIO_OPEN_OLD);
    }
    if (disk_qno_integrals_) {
        psio_->open(PSIF_DLPNO_QIA_QNO, PSIO_OPEN_NEW);
        psio_->open(PSIF_DLPNO_QAB_QNO, PSIO_OPEN_NEW);
    }

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

    if (disk_qno_integrals_) {
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
    set_scalar_variable("Q CORRECTION ENERGY",
                        e_ccsdtq_total - variables_["CCSDT TOTAL ENERGY"]);

    print_results();

    return e_ccsdtq_total;
}

void DLPNOCCSDTQ::print_results() {
    const int naocc = i_j_to_ij_.size();
    double t1_diagnostic = 0.0;
#pragma omp parallel for reduction(+ : t1_diagnostic)
    for (int i = 0; i < naocc; ++i) {
        t1_diagnostic += T_ia_[i]->vector_dot(T_ia_[i]);
    }
    t1_diagnostic = std::sqrt(t1_diagnostic / (2.0 * naocc));
    outfile->Printf("\n  T1 Diagnostic: %8.8f \n", t1_diagnostic);
    set_scalar_variable("CC T1 DIAGNOSTIC", t1_diagnostic);

    const double total_correlation =
        e_lccsdtq_ + de_lccsdt_q_screened_ + de_lccsd_t_screened_ + de_weak_ +
        de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
    const double total_energy = variables_["SCF TOTAL ENERGY"] + total_correlation;
    const double ccsdtq_minus_ccsdt = total_energy - variables_["CCSDT TOTAL ENERGY"];
    const double ccsdtq_minus_ccsdt_q = total_energy - variables_["CCSDT(Q) TOTAL ENERGY"];

    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSDTQ Correlation Energy: %16.12f \n", total_correlation);
    outfile->Printf("    LCCSDTQ Correlation Energy:           %16.12f \n", e_lccsdtq_);
    outfile->Printf("    Screened Triplets Contribution:       %16.12f \n", de_lccsd_t_screened_);
    outfile->Printf("    Screened Quadruplets Contribution:    %16.12f \n", de_lccsdt_q_screened_);
    outfile->Printf("    CCSDTQ - CCSDT Energy:                %16.12f \n", ccsdtq_minus_ccsdt);
    outfile->Printf("    CCSDTQ - CCSDT(Q) Energy:             %16.12f \n", ccsdtq_minus_ccsdt_q);
    outfile->Printf("\n\n  @Total DLPNO-CCSDTQ Energy: %16.12f \n", total_energy);
}

} // namespace dlpno
} // namespace psi
