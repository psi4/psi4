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

#include "streamed.h"

#include <algorithm>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/onebody.h"
#include "psi4/lib3index/dftensor.h"
#include <Einsums/TensorAlgebra.hpp>

namespace psi {
namespace f12 {

StreamedMP2F12::StreamedMP2F12(SharedWavefunction ref_wfn, Options& options)
    : MP2F12(ref_wfn, options), auxiliary_block_size_(options.get_int("F12_AUX_BLOCK_SIZE")) {}

void StreamedMP2F12::form_df_fock_streamed(einsums::Tensor<double, 2> *f,
                                        einsums::Tensor<double, 2> *k,
                                        einsums::Tensor<double, 2> *fk,
                                        const einsums::Tensor<double, 2> *metric_inverse,
                                        einsums::Tensor<double, 3> *raw_g_occ_ri) {
    using namespace einsums;
    using namespace einsums::tensor_algebra;
    using namespace einsums::index;

    MP2F12::form_oeints(f);

    if (raw_g_occ_ri == nullptr ||
        raw_g_occ_ri->dim(0) != static_cast<size_t>(naux_) ||
        raw_g_occ_ri->dim(1) != static_cast<size_t>(nocc_) ||
        raw_g_occ_ri->dim(2) != static_cast<size_t>(nri_)) {
        throw PSIEXCEPTION("Fock/pair shared G data has invalid dimensions");
    }
    form_oper_ints_pack({"G"}, {raw_g_occ_ri});

    Tensor<double, 1> j_trace{"B", naux_};
    {
        Tensor Id = create_identity_tensor("I", nocc_, nocc_);
        Tensor J_Oper = (*raw_g_occ_ri)(Range{0, naux_}, Range{0, nocc_}, Range{0, nocc_});
        einsum(Indices{B}, &j_trace, Indices{B, i, j}, J_Oper, Indices{i, j}, Id);
    }

    if (metric_inverse == nullptr) {
        throw PSIEXCEPTION("streamed Fock requires the shared DF metric inverse");
    }

    // For J-only blocks, reassociate the exact contraction
    //   sum_A j_trace(A) sum_B J^-1(A,B) BPQ(B,P,Q)
    // into one auxiliary weight.  This permits contraction in AO space before
    // the MO transform and removes the naux*nmo1*nmo2 metric materialization.
    Tensor<double, 1> j_weight{"direct-J auxiliary weight", naux_};
    einsum(Indices{B}, &j_weight, Indices{A, B}, *metric_inverse,
           Indices{A}, j_trace);

    outfile->Printf("     Forming J/K and reusable raw G context in blocked OBS/CABS blocks\n");
    const int spaces[4][2] = {{0, 0}, {0, 1}, {1, 0}, {1, 1}};
    for (const auto &space : spaces) {
        const int o1 = space[0];
        const int o2 = space[1];
        const auto nbf1 = bs_[o1].basisset()->nbf();
        const auto nbf2 = bs_[o2].basisset()->nbf();
        const auto nmo1 = o1 ? ncabs_ : nobs_;
        const auto nmo2 = o2 ? ncabs_ : nobs_;
        const auto off1 = o1 ? nobs_ : 0;
        const auto off2 = o2 ? nobs_ : 0;

        if (o2 == 1) {
            // CABS-ket blocks never contribute to K.  Generate one AO shell
            // pair at a time and reduce its auxiliary dimension immediately;
            // the full naux*nbf1*nbf2 AO tensor never exists.
            Tensor<double, 2> j_ao{"integral-direct-J AO block", nbf1, nbf2};
            three_index_direct_j_ao_computer(j_weight, &j_ao,
                                             bs_[o1].basisset(), bs_[o2].basisset());

            Tensor<double, 2> C1{"direct-J C1", nbf1, nmo1};
            Tensor<double, 2> C2{"direct-J C2", nbf2, nmo2};
            convert_C(&C1, bs_[o1], nbf1, nmo1, false);
            convert_C(&C2, bs_[o2], nbf2, nmo2, false);
            Tensor<double, 2> j_pQ{"direct-J pQ", nbf1, nmo2};
            Tensor<double, 2> j_PQ{"direct-J PQ", nmo1, nmo2};
            linear_algebra::gemm<false, false>(1.0, j_ao, C2, 0.0, &j_pQ);
            linear_algebra::gemm<true, false>(1.0, C1, j_pQ, 0.0, &j_PQ);

            for (int pidx = 0; pidx < nmo1; ++pidx) {
                for (int qidx = 0; qidx < nmo2; ++qidx) {
                    (*f)(off1 + pidx, off2 + qidx) += 2.0 * j_PQ(pidx, qidx);
                }
            }
            continue;
        }

        auto BPQ = std::make_unique<Tensor<double, 3>>("Streamed metric MO", naux_, nmo1, nmo2);
        Tensor<double, 2> C1{"Fock C1", nbf1, nmo1};
        Tensor<double, 2> C2{"Fock C2", nbf2, nmo2};
        convert_C(&C1, bs_[o1], nbf1, nmo1, false);
        convert_C(&C2, bs_[o2], nbf2, nmo2, false);
        three_index_mo_aux_blocked_pack({"G"}, {BPQ.get()}, bs_[o1].basisset(),
                                   bs_[o2].basisset(), C1, C2);

        Tensor<double, 3> metric_block{"Streamed orthogonalized metric", naux_, nmo1, nmo2};
        einsum(Indices{A, P, Q}, &metric_block, Indices{A, B}, *metric_inverse,
               Indices{B, P, Q}, *BPQ);
        BPQ.reset();

        Tensor f_block = (*f)(Range{off1, off1 + nmo1}, Range{off2, off2 + nmo2});
        einsum(1.0, Indices{P, Q}, &f_block, 2.0, Indices{B, P, Q}, metric_block,
               Indices{B}, j_trace);
        for (int pidx = 0; pidx < nmo1; pidx++) {
            for (int qidx = 0; qidx < nmo2; qidx++) {
                (*f)(off1 + pidx, off2 + qidx) = f_block(pidx, qidx);
            }
        }

        if (o2 == 0) {
            Tensor<double, 3> k_metric{"Streamed K metric", naux_, nocc_, nmo1};
            Tensor metric_occ = metric_block(All, All, Range{0, nocc_});
            permute(Indices{B, i, P}, &k_metric, Indices{B, P, i}, metric_occ);
            Tensor<double, 2> k_rows{"Streamed K rows", nmo1, nri_};
            Tensor K_Oper = (*raw_g_occ_ri)(Range{0, naux_}, Range{0, nocc_}, Range{0, nri_});
            einsum(Indices{P, Q}, &k_rows, Indices{B, i, P}, k_metric,
                   Indices{B, i, Q}, K_Oper);
            for (int pidx = 0; pidx < nmo1; pidx++) {
                for (int qidx = 0; qidx < nri_; qidx++) {
                    (*k)(off1 + pidx, qidx) = k_rows(pidx, qidx);
                }
            }
        }
    }

    for (int pidx = 0; pidx < nri_; pidx++) {
        for (int qidx = 0; qidx < nri_; qidx++) {
            const double f_before_exchange = (*f)(pidx, qidx);
            (*fk)(pidx, qidx) = f_before_exchange;
            (*f)(pidx, qidx) = f_before_exchange - (*k)(pidx, qidx);
        }
    }
}

double StreamedMP2F12::compute_energy() {
    using namespace einsums;
    using namespace einsums::linear_algebra;
    using namespace einsums::tensor_algebra;
    using namespace einsums::index;

    timer_on("MP2-F12 Compute Energy");
    einsums::profile::initialize();
    print_header();

    timer_on("OBS and CABS");
    form_basissets();
    timer_off("OBS and CABS");
    if (ncabs_ == 0) {
        throw PSIEXCEPTION("F12_SUBTYPE=STREAMED requires a nonempty CABS space");
    }

    outfile->Printf("\n ===> Streamed DF-MP2-F12 <===\n");
    outfile->Printf("     Auxiliary shell blocks and reusable occupied-pair workspaces\n\n");

    Tensor<double, 2> f{"Streamed Fock Matrix", nri_, nri_};
    Tensor<double, 2> k{"Streamed Exchange Matrix", nri_, nri_};
    Tensor<double, 2> fk{"Streamed Fock-before-exchange Matrix", nri_, nri_};
    Tensor<double, 2> metric_inverse{"shared JinvAB", naux_, naux_};
    auto raw_g_occ_ri = std::make_unique<Tensor<double, 3>>(
        "shared raw G occupied-RI", naux_, nocc_, nri_);
    for (int pidx = 0; pidx < nri_; ++pidx) {
        for (int qidx = 0; qidx < nri_; ++qidx) {
            f(pidx, qidx) = 0.0;
            k(pidx, qidx) = 0.0;
            fk(pidx, qidx) = 0.0;
        }
    }

    timer_on("shared DF metric inverse");
    form_metric_inverse(&metric_inverse);
    timer_off("shared DF metric inverse");

    timer_on("blocked-AO streamed Fock Matrix and shared G");
    form_df_fock_streamed(&f, &k, &fk, &metric_inverse, raw_g_occ_ri.get());
    timer_off("blocked-AO streamed Fock Matrix and shared G");

    timer_on("shared metric/operator contexts");
    // The Fock producer already generated the complete raw G occupied-by-RI
    // context.  Its active occupied rows are the G operator needed by the pair contractions;
    // apply the already-shared J^{-1} once to obtain the robust-DF metric.
    // This removes both duplicate CABS-RI insertions and their AO-to-MO work.
    Tensor<double, 3> oper_g{"shared G pair operator", naux_, nact_, nri_};
    for (int aux = 0; aux < naux_; ++aux) {
        for (int active = 0; active < nact_; ++active) {
            for (int orbital = 0; orbital < nri_; ++orbital) {
                oper_g(aux, active, orbital) =
                    (*raw_g_occ_ri)(aux, active + nfrzn_, orbital);
            }
        }
    }
    raw_g_occ_ri.reset();
    Tensor<double, 3> metric{"metric from shared G", naux_, nact_, nri_};
    einsum(Indices{A, P, Q}, &metric, Indices{A, B}, metric_inverse,
           Indices{B, P, Q}, oper_g);

    Tensor<double, 3> oper_f{"F operator", naux_, nact_, nri_};
    Tensor<double, 3> oper_f2{"F2 operator", naux_, nact_, nri_};
    Tensor<double, 3> oper_fg{"FG operator", naux_, nact_, nact_};
    Tensor<double, 3> oper_uf{"Uf operator", naux_, nact_, nact_};
    form_oper_ints_pack(
        std::vector<std::string>{"F", "F2"},
        std::vector<Tensor<double, 3>*>{&oper_f, &oper_f2});
    form_oper_ints_pack(
        std::vector<std::string>{"FG", "Uf"},
        std::vector<Tensor<double, 3>*>{&oper_fg, &oper_uf});

    Tensor<double, 2> aux_f{"F auxiliary operator", naux_, naux_};
    Tensor<double, 2> aux_f2{"F2 auxiliary operator", naux_, naux_};
    Tensor<double, 2> aux_fg{"FG auxiliary operator", naux_, naux_};
    Tensor<double, 2> aux_uf{"Uf auxiliary operator", naux_, naux_};
    form_oper_ints_pack(
        std::vector<std::string>{"F", "F2", "FG", "Uf"},
        std::vector<Tensor<double, 2>*>{&aux_f, &aux_f2, &aux_fg, &aux_uf});
    timer_off("shared metric/operator contexts");

    const std::vector<char> order_g = {'o', 'O', 'o', 'O', 'o', 'O', 'o', 'C'};
    const std::vector<char> order_f = {'o', 'O', 'o', 'O', 'o', 'O', 'o', 'C',
                                       'o', 'C', 'o', 'O', 'o', 'C', 'o', 'C'};
    const std::vector<char> order_f2 = {'o', 'o', 'o', 'O', 'o', 'o', 'o', 'C'};
    const std::vector<char> order_oo = {'o', 'o', 'o', 'o'};

    auto dot2 = [](const auto& lhs, const auto& rhs) {
        double value = 0.0;
        for (size_t row = 0; row < lhs.dim(0); ++row) {
            for (size_t col = 0; col < lhs.dim(1); ++col) value += lhs(row, col) * rhs(row, col);
        }
        return value;
    };
    auto sandwich_gemm = [&dot2](auto& lhs, auto& middle, auto& rhs,
                                 Tensor<double, 2>& product) {
        linear_algebra::gemm<false, false>(1.0, lhs, middle, 0.0, &product);
        return dot2(product, rhs);
    };
    auto right_transposed_gemm = [&dot2](auto& lhs, auto& rhs, auto& target,
                                         Tensor<double, 2>& product) {
        linear_algebra::gemm<false, true>(1.0, lhs, rhs, 0.0, &product);
        return dot2(product, target);
    };

    double E_f12_s = 0.0;
    double E_f12_t = 0.0;

    std::vector<std::pair<int, int>> occupied_pairs;
    occupied_pairs.reserve(nact_ * (nact_ + 1) / 2);
    for (int i = 0; i < nact_; ++i) {
        for (int j = i; j < nact_; ++j) occupied_pairs.emplace_back(i, j);
    }
    std::vector<PairResult> pair_results(occupied_pairs.size());
    const int pair_threads = std::min(nthreads_, static_cast<int>(occupied_pairs.size()));

    const double workspace_doubles =
        2.0 * nobs_ * nri_ + 2.0 * nri_ * nri_ + 2.0 * nact_ * nri_ +
        2.0 * nact_ * nact_ + 3.0 * nvir_ * nvir_ + nri_ * nri_ +
        nocc_ * nri_ + ncabs_ * nocc_ + nvir_ * nobs_ + 3.0 * naux_ * nri_;
    const double workspace_bytes = workspace_doubles * sizeof(double);
    const double total_workspace_bytes = workspace_bytes * pair_threads;
    if (total_workspace_bytes > Process::environment.get_memory()) {
        throw PSIEXCEPTION("STREAMED pair workspaces exceed the available memory; reduce the thread count");
    }

    std::vector<std::unique_ptr<PairWorkspace>> pair_workspaces;
    pair_workspaces.reserve(pair_threads);
    for (int workspace_id = 0; workspace_id < pair_threads; ++workspace_id) {
        pair_workspaces.push_back(std::make_unique<PairWorkspace>(
            workspace_id, nobs_, nri_, nocc_, nact_, nvir_, ncabs_, naux_));
    }

    outfile->Printf("     Pair workers: %d; workspace: %.3f MiB per worker\n",
                    pair_threads, workspace_bytes / (1024.0 * 1024.0));

    outfile->Printf("  %1s   %1s  |     %14s     %14s     %12s \n", "i", "j", "E_F12(Singlet)",
                    "E_F12(Triplet)", "E_F12");
    outfile->Printf(" ----------------------------------------------------------------------\n");
    timer_on("Streamed pair contraction");
#pragma omp parallel for schedule(static) num_threads(pair_threads)
    for (size_t pair_index = 0; pair_index < occupied_pairs.size(); ++pair_index) {
#ifdef _OPENMP
        const int workspace_id = omp_get_thread_num();
#else
        const int workspace_id = 0;
#endif
        auto& workspace = *pair_workspaces[workspace_id];
        const int i = occupied_pairs[pair_index].first;
        const int j = occupied_pairs[pair_index].second;
        auto& pair_result = pair_results[pair_index];
        // Each worker reuses its operator buffers for successive pairs.
        form_df_pair_ints_into("G", &metric, &oper_g, nullptr, i, j,
                               nobs_, nri_, order_g, &workspace.Gij, &workspace.df_scratch);
        form_df_pair_ints_into("F", &metric, &oper_f, &aux_f, i, j,
                               nri_, nri_, order_f, &workspace.Fij, &workspace.df_scratch);
        if (i != j) {
            // Electron interchange symmetry gives F_ji(q,s)=F_ij(s,q).
            for (int qidx = 0; qidx < nri_; ++qidx) {
                for (int sidx = 0; sidx < nri_; ++sidx) {
                    workspace.Fji(qidx, sidx) = workspace.Fij(sidx, qidx);
                }
            }
        }
        form_df_pair_ints_into("F2", &metric, &oper_f2, &aux_f2, i, j,
                               nact_, nri_, order_f2, &workspace.F2ij, &workspace.df_scratch);
        if (i != j) {
            form_df_pair_ints_into("G", &metric, &oper_g, nullptr, j, i,
                                   nobs_, nri_, order_g, &workspace.Gji, &workspace.df_scratch);
        }
        if (i != j) {
            form_df_pair_ints_into("F2", &metric, &oper_f2, &aux_f2, j, i,
                                   nact_, nri_, order_f2, &workspace.F2ji, &workspace.df_scratch);
        }
        form_df_pair_ints_into("FG", &metric, &oper_fg, &aux_fg, i, j,
                               nact_, nact_, order_oo, &workspace.FGij, &workspace.df_scratch);
        form_df_pair_ints_into("Uf", &metric, &oper_uf, &aux_uf, i, j,
                               nact_, nact_, order_oo, &workspace.Ufij, &workspace.df_scratch);

        auto& Gij = workspace.Gij;
        auto& Fij = workspace.Fij;
        auto& F2ij = workspace.F2ij;
        auto& Gji = (i == j) ? workspace.Gij : workspace.Gji;
        auto& Fji = (i == j) ? workspace.Fij : workspace.Fji;
        auto& F2ji =
            (i == j) ? workspace.F2ij : workspace.F2ji;
        auto& FGij = workspace.FGij;
        auto& Ufij = workspace.Ufij;

        auto Gij_oc = Gij(Range{0, nocc_}, Range{nobs_, nri_});
        auto Gji_oc = Gji(Range{0, nocc_}, Range{nobs_, nri_});
        auto Fij_oc = Fij(Range{0, nocc_}, Range{nobs_, nri_});
        auto Fji_oc = Fji(Range{0, nocc_}, Range{nobs_, nri_});
        auto Gij_pq = Gij(Range{0, nobs_}, Range{0, nobs_});
        auto Gji_pq = Gji(Range{0, nobs_}, Range{0, nobs_});
        auto Fij_pq = Fij(Range{0, nobs_}, Range{0, nobs_});
        auto Fji_pq = Fji(Range{0, nobs_}, Range{0, nobs_});

        const double v_ij = FGij(i, j) - dot2(Gij_oc, Fij_oc) - dot2(Gji_oc, Fji_oc) -
                            dot2(Fij_pq, Gij_pq);
        const double v_ji = (i == j)
                                ? v_ij
                                : FGij(j, i) - dot2(Gij_oc, Fji_oc) - dot2(Gji_oc, Fij_oc) -
                                      dot2(Fij_pq, Gji_pq);
        // F2's third occupied index is active-local, while its final RI
        // index retains the full occupied offset.
        const double x_ij = F2ij(i, j + nfrzn_) - dot2(Fij_oc, Fij_oc) - dot2(Fji_oc, Fji_oc) -
                            dot2(Fij_pq, Fij_pq);
        const double x_ji = (i == j)
                                ? x_ij
                                : F2ij(j, i + nfrzn_) - 2.0 * dot2(Fij_oc, Fji_oc) -
                                      dot2(Fij_pq, Fji_pq);

        auto Fij_vc = Fij(Range{nocc_, nobs_}, Range{nobs_, nri_});
        auto Fji_vc = Fji(Range{nocc_, nobs_}, Range{nobs_, nri_});
        auto f_vc = f(Range{nocc_, nobs_}, Range{nobs_, nri_});
        auto& Cij = workspace.Cij;
        auto& Cji = workspace.Cji;
        linear_algebra::gemm<false, true>(1.0, Fij_vc, f_vc, 0.0, &Cij);
        linear_algebra::gemm<false, true>(1.0, f_vc, Fji_vc, 1.0, &Cij);
        if (i == j) {
            Cji = Cij;
        } else {
            linear_algebra::gemm<false, true>(1.0, Fji_vc, f_vc, 0.0, &Cji);
            linear_algebra::gemm<false, true>(1.0, f_vc, Fij_vc, 1.0, &Cji);
        }

        auto& Dij = workspace.Dij;
        for (int aidx = 0; aidx < nvir_; ++aidx) {
            for (int bidx = 0; bidx < nvir_; ++bidx) {
                Dij(aidx, bidx) = 1.0 /
                    (f(nocc_ + aidx, nocc_ + aidx) + f(nocc_ + bidx, nocc_ + bidx) -
                     f(nfrzn_ + i, nfrzn_ + i) - f(nfrzn_ + j, nfrzn_ + j));
            }
        }

        double b_ij = Ufij(i, j);
        double b_ji = (i == j) ? b_ij : Ufij(j, i);
        for (int A = 0; A < nri_; ++A) {
            b_ij += fk(nfrzn_ + i, A) * F2ji(j, A) + F2ij(i, A) * fk(nfrzn_ + j, A);
            if (i != j) {
                b_ji += fk(nfrzn_ + i, A) * F2ij(j, A) + F2ji(i, A) * fk(nfrzn_ + j, A);
            }
        }

        b_ij -= sandwich_gemm(Fij, k, Fij, workspace.product_full) +
                sandwich_gemm(Fji, k, Fji, workspace.product_full);
        if (i != j) {
            b_ji -= sandwich_gemm(Fij, k, Fji, workspace.product_full) +
                    sandwich_gemm(Fji, k, Fij, workspace.product_full);
        }

        auto Fij_o1 = Fij(Range{0, nocc_}, All);
        auto Fji_o1 = Fji(Range{0, nocc_}, All);
        b_ij -= sandwich_gemm(Fij_o1, f, Fij_o1, workspace.product_occ_all) +
                sandwich_gemm(Fji_o1, f, Fji_o1, workspace.product_occ_all);
        if (i != j) {
            b_ji -= sandwich_gemm(Fij_o1, f, Fji_o1, workspace.product_occ_all) +
                    sandwich_gemm(Fji_o1, f, Fij_o1, workspace.product_occ_all);
        }

        auto Fij_co = Fij(Range{nobs_, nri_}, Range{0, nocc_});
        auto Fji_co = Fji(Range{nobs_, nri_}, Range{0, nocc_});
        auto f_oo = f(Range{0, nocc_}, Range{0, nocc_});
        double b57_ij = sandwich_gemm(Fij_co, f_oo, Fij_co, workspace.product_cabs_occ) +
                        sandwich_gemm(Fji_co, f_oo, Fji_co, workspace.product_cabs_occ);
        double b57_ji = 0.0;
        if (i != j) {
            b57_ji = sandwich_gemm(Fij_co, f_oo, Fji_co, workspace.product_cabs_occ) +
                     sandwich_gemm(Fji_co, f_oo, Fij_co, workspace.product_cabs_occ);
        }
        auto Fij_c1 = Fij(Range{nobs_, nri_}, All);
        auto Fji_c1 = Fji(Range{nobs_, nri_}, All);
        auto f_o1 = f(Range{0, nocc_}, All);
        b57_ij -= 2.0 * right_transposed_gemm(Fij_c1, f_o1, Fij_co, workspace.product_cabs_occ) +
                  2.0 * right_transposed_gemm(Fji_c1, f_o1, Fji_co, workspace.product_cabs_occ);
        if (i != j) {
            b57_ji -= 2.0 * right_transposed_gemm(Fij_c1, f_o1, Fji_co, workspace.product_cabs_occ) +
                      2.0 * right_transposed_gemm(Fji_c1, f_o1, Fij_co, workspace.product_cabs_occ);
        }
        b_ij += b57_ij;
        if (i != j) b_ji += b57_ji;

        auto Fij_vq = Fij(Range{nocc_, nobs_}, Range{0, nobs_});
        auto Fji_vq = Fji(Range{nocc_, nobs_}, Range{0, nobs_});
        auto f_pq = f(Range{0, nobs_}, Range{0, nobs_});
        double b68_ij = sandwich_gemm(Fij_vq, f_pq, Fij_vq, workspace.product_vir_obs) +
                        sandwich_gemm(Fji_vq, f_pq, Fji_vq, workspace.product_vir_obs);
        double b68_ji = 0.0;
        if (i != j) {
            b68_ji = sandwich_gemm(Fij_vq, f_pq, Fji_vq, workspace.product_vir_obs) +
                     sandwich_gemm(Fji_vq, f_pq, Fij_vq, workspace.product_vir_obs);
        }
        auto f_pc = f(Range{0, nobs_}, Range{nobs_, nri_});
        b68_ij += 2.0 * right_transposed_gemm(Fij_vc, f_pc, Fij_vq, workspace.product_vir_obs) +
                  2.0 * right_transposed_gemm(Fji_vc, f_pc, Fji_vq, workspace.product_vir_obs);
        if (i != j) {
            b68_ji += 2.0 * right_transposed_gemm(Fij_vc, f_pc, Fji_vq, workspace.product_vir_obs) +
                      2.0 * right_transposed_gemm(Fji_vc, f_pc, Fij_vq, workspace.product_vir_obs);
        }
        b_ij -= b68_ij;
        if (i != j) b_ji -= b68_ji;
        if (i == j) b_ji = b_ij;

        const double eps_sum = f(nfrzn_ + i, nfrzn_ + i) + f(nfrzn_ + j, nfrzn_ + j);
        const double be_ij = b_ij - eps_sum * x_ij;
        const double be_ji = b_ji - eps_sum * x_ji;

        auto Gij_vv = Gij(Range{nocc_, nobs_}, Range{nocc_, nobs_});
        double cdg_ij = 0.0, cdg_ji = 0.0, cdc_ij = 0.0, cdc_ji = 0.0;
        for (int aidx = 0; aidx < nvir_; ++aidx) {
            for (int bidx = 0; bidx < nvir_; ++bidx) {
                const double gd = Gij_vv(aidx, bidx) * Dij(aidx, bidx);
                const double cd = Cij(aidx, bidx) * Dij(aidx, bidx);
                cdg_ij += Cij(aidx, bidx) * gd;
                cdc_ij += Cij(aidx, bidx) * cd;
                if (i != j) {
                    cdg_ji += Cji(aidx, bidx) * gd;
                    cdc_ji += Cji(aidx, bidx) * cd;
                }
            }
        }
        if (i == j) {
            cdg_ji = cdg_ij;
            cdc_ji = cdc_ij;
        }

        const double vt_ij = v_ij - cdg_ij;
        const double vt_ji = v_ji - cdg_ji;
        const double bt_ij = be_ij - cdc_ij;
        const double bt_ji = be_ji - cdc_ji;
        const int kd = (i == j) ? 1 : 2;
        const double t_plus = t_(i, j, i, j) + t_(i, j, j, i);
        const double vt_s = 0.25 * t_plus * kd * (vt_ij + vt_ji);
        const double bt_s = 0.125 * t_plus * kd * (bt_ij + bt_ji) * t_plus * kd;
        const double E_s = kd * (2.0 * vt_s + bt_s);
        pair_result.energy_s = E_s;

        double E_t = 0.0;
        if (i != j) {
            const double t_minus = t_(i, j, i, j) - t_(i, j, j, i);
            const double vt_t = 0.25 * t_minus * kd * (vt_ij - vt_ji);
            const double bt_t = 0.125 * t_minus * kd * (bt_ij - bt_ji) * t_minus * kd;
            E_t = 3.0 * kd * (2.0 * vt_t + bt_t);
            pair_result.energy_t = E_t;
        }
    }
    timer_off("Streamed pair contraction");

    // Sum in occupied-pair order, independently of thread completion order.
    for (size_t pair_index = 0; pair_index < occupied_pairs.size(); ++pair_index) {
        const auto& pair_result = pair_results[pair_index];
        const int i = occupied_pairs[pair_index].first;
        const int j = occupied_pairs[pair_index].second;
        E_f12_s += pair_result.energy_s;
        E_f12_t += pair_result.energy_t;
        outfile->Printf("%3d %3d  |   %16.12f   %16.12f     %16.12f \n", i + nfrzn_ + 1, j + nfrzn_ + 1,
                        pair_result.energy_s, pair_result.energy_t,
                        pair_result.energy_s + pair_result.energy_t);
    }

    set_scalar_variable("MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY",
                        E_f12_s + scalar_variable("MP2 OPPOSITE-SPIN CORRELATION ENERGY"));
    set_scalar_variable("MP2-F12 SAME-SPIN CORRELATION ENERGY",
                        E_f12_t + scalar_variable("MP2 SAME-SPIN CORRELATION ENERGY"));
    E_f12_ = E_f12_s + E_f12_t;

    if (singles_) {
        timer_on("CABS Singles Correction");
        MP2F12::form_cabs_singles(&f);
        timer_off("CABS Singles Correction");
    }
    set_scalar_variable("F12 CABS CORRECTION ENERGY", E_singles_);
    print_results();
    set_energy(E_mp2f12_);
    if (print_ > 1) einsums::profile::report("timer_mp2f12_streamed.dat", false);
    einsums::profile::finalize();
    timer_off("MP2-F12 Compute Energy");
    return E_mp2f12_;
}

}  // namespace f12
}  // namespace psi
