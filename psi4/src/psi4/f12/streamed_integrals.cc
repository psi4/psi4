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
#include <cstring>
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

void StreamedMP2F12::three_index_direct_j_ao_computer(
    const einsums::Tensor<double, 1>& j_weight,
    einsums::Tensor<double, 2>* j_ao,
    std::shared_ptr<BasisSet> bs1,
    std::shared_ptr<BasisSet> bs2) {
    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    std::shared_ptr<IntegralFactory> intf(new IntegralFactory(DFBS_, zero, bs1, bs2));
    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    ints.push_back(std::shared_ptr<TwoBodyAOInt>(intf->eri()));
    for (size_t thread = 1; thread < nthreads_; ++thread) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
    }

    const bool equivalent_bases = bs1 == bs2;

    // Assign one AO shell pair to exactly one thread.  Auxiliary shells are
    // consumed serially inside that task, so no atomics or thread-private
    // nbf*nbf matrices are required.
#pragma omp parallel for collapse(2) schedule(guided) num_threads(nthreads_)
    for (size_t P = 0; P < bs1->nshell(); ++P) {
        for (size_t Q = 0; Q < bs2->nshell(); ++Q) {
            if (equivalent_bases && Q < P) continue;

            size_t rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            const auto numP = bs1->shell(P).nfunction();
            const auto numQ = bs2->shell(Q).nfunction();
            const auto index_P = bs1->shell(P).function_index();
            const auto index_Q = bs2->shell(Q).function_index();

            for (size_t p = 0; p < numP; ++p) {
                for (size_t q = 0; q < numQ; ++q) {
                    (*j_ao)(index_P + p, index_Q + q) = 0.0;
                }
            }

            for (size_t B = 0; B < DFBS_->nshell(); ++B) {
                const auto numB = DFBS_->shell(B).nfunction();
                const auto index_B = DFBS_->shell(B).function_index();
                ints[rank]->compute_shell(B, 0, P, Q);
                const auto* buffer = ints[rank]->buffers()[0];

                for (size_t b = 0, offset = 0; b < numB; ++b) {
                    const double weight = j_weight(index_B + b);
                    for (size_t p = 0; p < numP; ++p) {
                        for (size_t q = 0; q < numQ; ++q, ++offset) {
                            (*j_ao)(index_P + p, index_Q + q) +=
                                weight * buffer[offset];
                        }
                    }
                }
            }

            if (equivalent_bases && P != Q) {
                for (size_t p = 0; p < numP; ++p) {
                    for (size_t q = 0; q < numQ; ++q) {
                        (*j_ao)(index_Q + q, index_P + p) =
                            (*j_ao)(index_P + p, index_Q + q);
                    }
                }
            }
        }
    }
}

void StreamedMP2F12::three_index_mo_aux_blocked_pack(
    const std::vector<std::string>& int_types,
    const std::vector<einsums::Tensor<double, 3>*>& BPQ_outputs,
    std::shared_ptr<BasisSet> bs1,
    std::shared_ptr<BasisSet> bs2,
    const einsums::Tensor<double, 2>& C1,
    const einsums::Tensor<double, 2>& C2) {
    using namespace einsums;
    using namespace einsums::tensor_algebra;
    using namespace einsums::index;

    if (int_types.empty() || int_types.size() != BPQ_outputs.size()) {
        throw PSIEXCEPTION("operator pack types/outputs mismatch");
    }
    for (auto* output : BPQ_outputs) {
        if (output == nullptr || output->dim(0) != static_cast<size_t>(naux_) ||
            output->dim(1) != C1.dim(1) || output->dim(2) != C2.dim(1)) {
            throw PSIEXCEPTION("operator pack output dimension mismatch");
        }
    }

    const int target_functions = auxiliary_block_size_;

    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    std::shared_ptr<IntegralFactory> intf(new IntegralFactory(DFBS_, zero, bs1, bs2));
    std::vector<std::vector<std::shared_ptr<TwoBodyAOInt>>> evaluators(int_types.size());
    for (size_t op = 0; op < int_types.size(); ++op) {
        const auto& type = int_types[op];
        if (type == "F") {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12(cgtg_)));
        } else if (type == "FG") {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12g12(cgtg_)));
        } else if (type == "F2") {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_squared(cgtg_)));
        } else if (type == "Uf") {
            evaluators[op].push_back(
                std::shared_ptr<TwoBodyAOInt>(intf->f12_double_commutator(cgtg_)));
        } else {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->eri()));
        }
        for (size_t thread = 1; thread < nthreads_; ++thread) {
            evaluators[op].push_back(
                std::shared_ptr<TwoBodyAOInt>(evaluators[op][0]->clone()));
        }
    }

    const bool equivalent_bases = bs1 == bs2;
    size_t block_count = 0;
    size_t maximum_block_functions = 0;
    size_t shell_begin = 0;
    while (shell_begin < DFBS_->nshell()) {
        const size_t function_begin = DFBS_->shell(shell_begin).function_index();
        size_t shell_end = shell_begin;
        size_t block_functions = 0;
        while (shell_end < DFBS_->nshell()) {
            const size_t shell_functions = DFBS_->shell(shell_end).nfunction();
            if (block_functions > 0 &&
                block_functions + shell_functions > static_cast<size_t>(target_functions)) {
                break;
            }
            block_functions += shell_functions;
            ++shell_end;
        }

        std::vector<std::unique_ptr<Tensor<double, 3>>> ao_blocks;
        ao_blocks.reserve(int_types.size());
        for (const auto& type : int_types) {
            ao_blocks.push_back(std::make_unique<Tensor<double, 3>>(
                "" + type + " auxiliary AO block", block_functions,
                bs1->nbf(), bs2->nbf()));
        }

#pragma omp parallel for collapse(3) schedule(guided) num_threads(nthreads_)
        for (size_t Bshell = shell_begin; Bshell < shell_end; ++Bshell) {
            for (size_t P = 0; P < bs1->nshell(); ++P) {
                for (size_t Q = 0; Q < bs2->nshell(); ++Q) {
                    if (equivalent_bases && Q < P) continue;
                    size_t rank = 0;
#ifdef _OPENMP
                    rank = omp_get_thread_num();
#endif
                    const auto numB = DFBS_->shell(Bshell).nfunction();
                    const auto numP = bs1->shell(P).nfunction();
                    const auto numQ = bs2->shell(Q).nfunction();
                    const auto index_B = DFBS_->shell(Bshell).function_index();
                    const auto index_P = bs1->shell(P).function_index();
                    const auto index_Q = bs2->shell(Q).function_index();

                    for (size_t op = 0; op < int_types.size(); ++op) {
                        evaluators[op][rank]->compute_shell(Bshell, 0, P, Q);
                        const auto* buffer = evaluators[op][rank]->buffers()[0];
                        auto& Bpq = *ao_blocks[op];
                        for (size_t b = 0, offset = 0; b < numB; ++b) {
                            const size_t local_B = index_B + b - function_begin;
                            for (size_t p = 0; p < numP; ++p) {
                                const size_t function_P = index_P + p;
                                for (size_t q = 0; q < numQ; ++q, ++offset) {
                                    const size_t function_Q = index_Q + q;
                                    Bpq(local_B, function_P, function_Q) = buffer[offset];
                                    if (equivalent_bases && P != Q) {
                                        Bpq(local_B, function_Q, function_P) = buffer[offset];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (size_t op = 0; op < int_types.size(); ++op) {
            Tensor<double, 3> BpQ{"packed BpQ block", block_functions,
                                  bs1->nbf(), C2.dim(1)};
            einsum(Indices{B, p, Q}, &BpQ, Indices{B, p, q}, *ao_blocks[op],
                   Indices{q, Q}, C2);
            ao_blocks[op].reset();
            Tensor<double, 3> BPQ_block{"packed BPQ block", block_functions,
                                        C1.dim(1), C2.dim(1)};
            einsum(Indices{B, P, Q}, &BPQ_block, Indices{B, p, Q}, BpQ,
                   Indices{p, P}, C1);
            for (size_t B = 0; B < block_functions; ++B) {
                for (size_t P = 0; P < C1.dim(1); ++P) {
                    for (size_t Q = 0; Q < C2.dim(1); ++Q) {
                        (*BPQ_outputs[op])(function_begin + B, P, Q) =
                            BPQ_block(B, P, Q);
                    }
                }
            }
        }

        maximum_block_functions = std::max(maximum_block_functions, block_functions);
        ++block_count;
        shell_begin = shell_end;
    }

    outfile->Printf("     operator pack width %zu: %zu blocks, max %zu functions/block\n",
                    int_types.size(), block_count, maximum_block_functions);
}

void StreamedMP2F12::form_metric_inverse(
    einsums::Tensor<double, 2>* metric_inverse) {
    FittingMetric metric(DFBS_, true);
    metric.form_full_eig_inverse(1.0e-12);
    SharedMatrix Jinv = metric.get_metric();
    double** rows = Jinv->pointer();
    for (size_t A = 0; A < naux_; ++A) {
        std::memcpy(&(*metric_inverse)(A, 0), rows[A], naux_ * sizeof(double));
    }
}

void StreamedMP2F12::form_oper_ints_pack(
    const std::vector<std::string>& int_types,
    const std::vector<einsums::Tensor<double, 3>*>& DF_ERI_outputs) {
    using namespace einsums;

    if (int_types.empty() || int_types.size() != DF_ERI_outputs.size()) {
        throw PSIEXCEPTION("three-center operator pack types/outputs mismatch");
    }
    for (auto* output : DF_ERI_outputs) {
        if (output == nullptr) {
            throw PSIEXCEPTION("three-center operator pack requires non-null outputs");
        }
    }
    const auto dim1 = DF_ERI_outputs[0]->dim(1);
    const auto dim2 = DF_ERI_outputs[0]->dim(2);
    for (auto* output : DF_ERI_outputs) {
        if (output->dim(0) != static_cast<size_t>(naux_) ||
            output->dim(1) != dim1 || output->dim(2) != dim2) {
            throw PSIEXCEPTION("three-center operator pack requires identical output shapes");
        }
    }

    const std::vector<char> order = dim2 == static_cast<size_t>(nri_)
                                       ? std::vector<char>{'o', 'O', 'o', 'C'}
                                       : std::vector<char>{'o', 'o'};
    const bool use_offset = dim2 == static_cast<size_t>(nri_);
    const bool frzn_1 = dim1 == static_cast<size_t>(nact_);
    const bool frzn_2 = dim2 == static_cast<size_t>(nact_);

    for (int idx = 0; idx < static_cast<int>(order.size() / 2); ++idx) {
        const int i = idx * 2;
        const int o1 = order[i] == 'C' ? 1 : 0;
        const int o2 = order[i + 1] == 'C' ? 1 : 0;
        const auto nbf1 = bs_[o1].basisset()->nbf();
        const auto nbf2 = bs_[o2].basisset()->nbf();
        const auto nmo1 = dim1;
        const auto nmo2 = o2 ? static_cast<size_t>(ncabs_)
                              : (order[i + 1] == 'O' ? static_cast<size_t>(nobs_) : dim2);

        Tensor<double, 2> C1{"packed operator C1", nbf1, nmo1};
        Tensor<double, 2> C2{"packed operator C2", nbf2, nmo2};
        convert_C(&C1, bs_[o1], nbf1, nmo1, frzn_1);
        convert_C(&C2, bs_[o2], nbf2, nmo2, frzn_2);

        std::vector<std::unique_ptr<Tensor<double, 3>>> packed_blocks;
        std::vector<Tensor<double, 3>*> packed_block_ptrs;
        packed_blocks.reserve(int_types.size());
        packed_block_ptrs.reserve(int_types.size());
        for (const auto& type : int_types) {
            packed_blocks.push_back(std::make_unique<Tensor<double, 3>>(
                "" + type + " packed BPQ", naux_, nmo1, nmo2));
            packed_block_ptrs.push_back(packed_blocks.back().get());
        }
        three_index_mo_aux_blocked_pack(int_types, packed_block_ptrs,
                                        bs_[o1].basisset(), bs_[o2].basisset(), C1, C2);

        const auto off1 = o1 && use_offset ? nobs_ : 0;
        const auto off2 = o2 && use_offset ? nobs_ : 0;
        for (size_t op = 0; op < int_types.size(); ++op) {
            TensorView<double, 3> destination{
                *DF_ERI_outputs[op], Dim<3>{static_cast<size_t>(naux_), nmo1, nmo2},
                Offset<3>{0, static_cast<size_t>(off1), static_cast<size_t>(off2)}};
            destination = *packed_blocks[op];
        }
    }
}

void StreamedMP2F12::form_oper_ints_pack(
    const std::vector<std::string>& int_types,
    const std::vector<einsums::Tensor<double, 2>*>& DF_ERI_outputs) {
    if (int_types.empty() || int_types.size() != DF_ERI_outputs.size()) {
        throw PSIEXCEPTION("auxiliary operator pack types/outputs mismatch");
    }
    for (auto* output : DF_ERI_outputs) {
        if (output == nullptr || output->dim(0) != static_cast<size_t>(naux_) ||
            output->dim(1) != static_cast<size_t>(naux_)) {
            throw PSIEXCEPTION("auxiliary operator pack output dimension mismatch");
        }
    }

    std::shared_ptr<BasisSet> zero(BasisSet::zero_ao_basis_set());
    std::shared_ptr<IntegralFactory> intf(new IntegralFactory(DFBS_, zero, DFBS_, zero));
    std::vector<std::vector<std::shared_ptr<TwoBodyAOInt>>> evaluators(int_types.size());
    for (size_t op = 0; op < int_types.size(); ++op) {
        const auto& type = int_types[op];
        if (type == "F") {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12(cgtg_)));
        } else if (type == "FG") {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12g12(cgtg_)));
        } else if (type == "F2") {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->f12_squared(cgtg_)));
        } else if (type == "Uf") {
            evaluators[op].push_back(
                std::shared_ptr<TwoBodyAOInt>(intf->f12_double_commutator(cgtg_)));
        } else {
            evaluators[op].push_back(std::shared_ptr<TwoBodyAOInt>(intf->eri()));
        }
        for (size_t thread = 1; thread < nthreads_; ++thread) {
            evaluators[op].push_back(
                std::shared_ptr<TwoBodyAOInt>(evaluators[op][0]->clone()));
        }
    }

#pragma omp parallel for collapse(2) schedule(guided) num_threads(nthreads_)
    for (size_t A = 0; A < DFBS_->nshell(); ++A) {
        for (size_t B = 0; B < DFBS_->nshell(); ++B) {
            if (B > A) continue;

            size_t rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            const size_t numA = DFBS_->shell(A).nfunction();
            const size_t numB = DFBS_->shell(B).nfunction();
            const size_t index_A = DFBS_->shell(A).function_index();
            const size_t index_B = DFBS_->shell(B).function_index();
            for (size_t op = 0; op < int_types.size(); ++op) {
                evaluators[op][rank]->compute_shell(A, 0, B, 0);
                const auto* buffer = evaluators[op][rank]->buffers()[0];
                size_t offset = 0;
                for (size_t a = 0; a < numA; ++a) {
                    for (size_t b = 0; b < numB; ++b) {
                        const double value = buffer[offset++];
                        (*DF_ERI_outputs[op])(index_A + a, index_B + b) = value;
                        if (A != B) (*DF_ERI_outputs[op])(index_B + b, index_A + a) = value;
                    }
                }
            }
        }
    }
}

void StreamedMP2F12::form_df_pair_ints_into(
    const std::string& int_type, einsums::Tensor<double, 3>* metric,
    einsums::Tensor<double, 3>* oper, einsums::Tensor<double, 2>* aux_oper,
    int pidx, int ridx, int rows, int cols, const std::vector<char>& order,
    einsums::Tensor<double, 2>* result, PairDFScratch* scratch) {
    using namespace einsums;
    using namespace einsums::linear_algebra;
    using namespace einsums::tensor_algebra;
    using namespace einsums::index;

    if (result->dim(0) != static_cast<size_t>(rows) ||
        result->dim(1) != static_cast<size_t>(cols)) {
        throw PSIEXCEPTION("reusable pair-integral buffer has incompatible dimensions");
    }
    for (int qidx = 0; qidx < rows; qidx++) {
        for (int sidx = 0; sidx < cols; sidx++) (*result)(qidx, sidx) = 0.0;
    }

    const bool frz_ket1 = (nfrzn_ > 0) && (rows == nact_);
    const bool frz_ket2 = (nfrzn_ > 0) && (cols == nact_);
    for (int block_index = 0; block_index < static_cast<int>(order.size() / 4); block_index++) {
        const int offset = block_index * 4;
        const auto nmo2 = (order[offset + 1] == 'C')
                              ? ncabs_
                              : (order[offset + 1] == 'O') ? nobs_ : (frz_ket1 ? nact_ : nocc_);
        const auto nmo4 = (order[offset + 3] == 'C')
                              ? ncabs_
                              : (order[offset + 3] == 'O') ? nobs_ : (frz_ket2 ? nact_ : nocc_);
        const auto off2 = (order[offset + 1] == 'C') ? nobs_ : 0;
        const auto off4 = (order[offset + 3] == 'C') ? nobs_ : 0;

        const auto term1_off2 = frz_ket1 ? nfrzn_ : 0;
        const bool full_ri_operator = oper->dim(2) == static_cast<size_t>(nri_);
        const auto term1_off4 = (frz_ket2 && full_ri_operator) ? nfrzn_ : 0;
        auto left_metric = scratch->aux_left(All, Range{0, nmo2});
        auto right_oper = scratch->aux_right(All, Range{0, nmo4});
        for (int aux = 0; aux < naux_; aux++) {
            for (int qidx = 0; qidx < nmo2; qidx++) {
                left_metric(aux, qidx) =
                    (*metric)(aux, pidx, off2 + term1_off2 + qidx);
            }
            for (int sidx = 0; sidx < nmo4; sidx++) {
                right_oper(aux, sidx) =
                    (*oper)(aux, ridx, off4 + term1_off4 + sidx);
            }
        }

        // Each order entry addresses a disjoint result block, so the robust-DF
        // contraction can target the final pair buffer directly.
        auto block = (*result)(Range{off2, off2 + nmo2}, Range{off4, off4 + nmo4});
        linear_algebra::gemm<true, false>(1.0, left_metric, right_oper,
                                          0.0, &block);

        if (int_type != "G") {
            const auto term2_off2 = (frz_ket1 && full_ri_operator) ? nfrzn_ : 0;
            const auto term2_off4 = frz_ket2 ? nfrzn_ : 0;
            // Term 1 has consumed right_oper, so reuse that panel for
            // left_oper.  aux_other retains right_metric through term 3.
            auto left_oper = scratch->aux_right(All, Range{0, nmo2});
            auto right_metric = scratch->aux_other(All, Range{0, nmo4});
            for (int aux = 0; aux < naux_; aux++) {
                for (int qidx = 0; qidx < nmo2; qidx++) {
                    left_oper(aux, qidx) =
                        (*oper)(aux, pidx, off2 + term2_off2 + qidx);
                }
                for (int sidx = 0; sidx < nmo4; sidx++) {
                    right_metric(aux, sidx) =
                        (*metric)(aux, ridx, off4 + term2_off4 + sidx);
                }
            }
            linear_algebra::gemm<true, false>(1.0, left_oper, right_metric,
                                              1.0, &block);

            // left_oper is dead after term 2; overwrite the same panel with
            // aux_oper * right_metric and finish term 3 in-place.
            auto metric_oper_product = scratch->aux_right(All, Range{0, nmo4});
            linear_algebra::gemm<false, false>(1.0, *aux_oper, right_metric,
                                               0.0, &metric_oper_product);
            linear_algebra::gemm<true, false>(-1.0, left_metric,
                                              metric_oper_product, 1.0,
                                              &block);
        }
    }
}

}  // namespace f12
}  // namespace psi
