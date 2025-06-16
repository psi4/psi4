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
#include "dct_df_tensor.h"

#include "psi4/libpsi4util/process.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdpd/dpd.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dct {

DFTensor::DFTensor(const std::string& name, const int nQ, const Dimension& idx2, const Dimension& idx3) : Matrix(name, std::vector<int>(idx2.n(), nQ), setup(idx2, idx3))  {
    nQ_ = nQ;
    dim2_ = idx2;
    dim3_ = idx3;
}

DFTensor::DFTensor(const std::string& name, const DFTensor& tensor) {
    DFTensor(name, tensor.nQ(), tensor.idx2pi(), tensor.idx3pi());
}

const Dimension DFTensor::setup(const Dimension& idx2, const Dimension& idx3) {
    // Set up dimensions for b(Aux|PQ)
    if (idx2.n() != idx3.n())
        throw PSIEXCEPTION("DCT::DFTensor initialization fail. Inconsistent number of irreps for the two primary dimensions.");
    auto nirrep = idx2.n();
    Dimension LR(nirrep);
    for (int hL = 0; hL < nirrep; ++hL) {
        for (int hR = 0; hR < nirrep; ++hR) {
            LR[hL ^ hR] += idx2[hL] * idx3[hR];
        }
    }

    return LR;
}

// Convenience function for the simple case of a (Q|pq) pr qs -> (Q|qs).
DFTensor DFTensor::three_idx_primary_transform(const DFTensor& three_idx, const Matrix& left, const Matrix& right) {
    auto nQ = three_idx.rowdim(0);
    auto result = DFTensor("Three-Index Tensor", nQ, left.colspi(), right.colspi());
    result.three_idx_primary_transform_gemm(three_idx, left, right, 1.0, 0.0);

    return result;
}

// TODO: This should probably be migrated to/replaced with lib3index's DFHelper.
// However, we need symmetry, and lib3index currently doesn't support it. JPM 01/2021
void DFTensor::three_idx_primary_transform_gemm(const DFTensor& three_idx, const Matrix& left, const Matrix& right,
                                                 double alpha, double beta) {
    dct_timer_on("DCTSolver::Three-Index SO -> MO");

    if (three_idx.symmetry() || left.symmetry() || right.symmetry() || symmetry())
        throw PSIEXCEPTION("three_idx_primary_transform_gemm: Can only handle totally symmetric matrices.");

    if (three_idx.nirrep() != left.nirrep() || three_idx.nirrep() != right.nirrep() ||
        three_idx.nirrep() != nirrep()) {
        throw PSIEXCEPTION("three_idx_primary_transform_gemm: Number of irreps don't equal.");
    }

    if (three_idx.rowspi() != rowspi()) {
        throw PSIEXCEPTION(
            "three_idx_primary_transform_gemm: Tensor to transform and result must agree about number of number of "
            "aux. functions");
    }

    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    auto nQ = three_idx.rowdim(0);
    std::vector<int> offset_mo(three_idx.nirrep(), 0), offset_so(three_idx.nirrep(), 0);

    for (int h = 0; h < nirrep_; ++h) {
        auto three_idx_p = three_idx.pointer(h);
        auto result_p = pointer(h);
        for (int hL = 0; hL < nirrep_; ++hL) {
            const auto hR = h ^ hL;
            if (left.colspi(hL) > 0 && right.colspi(hR) > 0 && left.rowspi(hL) > 0 && right.rowspi(hR) > 0) {
                const auto leftP = left.pointer(hL);
                const auto rightP = right.pointer(hR);
                auto tmp = Matrix("Half-Transformed", nQ, left.rowspi(hL) * right.colspi(hR));
                auto tmpp = tmp.pointer();
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (int Q = 0; Q < nQ; ++Q) {
                    // First-half transformation
                    C_DGEMM('N', 'N', left.rowspi(hL), right.colspi(hR), right.rowspi(hR), 1.0,
                            three_idx_p[Q] + offset_so[h], right.rowspi(hR), rightP[0], right.colspi(hR), 0.0, tmpp[Q],
                            right.colspi(hR));
                    // Second-half transformation
                    C_DGEMM('T', 'N', left.colspi(hL), right.colspi(hR), left.rowspi(hL), alpha, leftP[0],
                            left.colspi(hL), tmpp[Q], right.colspi(hR), beta, result_p[Q] + offset_mo[h],
                            right.colspi(hR));
                }
            }
            offset_so[h] += left.rowspi(hL) * right.rowspi(hR);
            offset_mo[h] += left.colspi(hL) * right.colspi(hR);
        }
        if (offset_so[h] != three_idx.colspi(h)) throw PSIEXCEPTION("three_idx_primary_transform: Dimension mismatch");
    }

    dct_timer_off("DCTSolver::Three-Index SO -> MO");
}

void DFTensor::contract343(const DFTensor& b, dpdbuf4& G, bool transpose, double alpha, double beta) {
    if (b.rowspi() != rowspi_) {
        throw PSIEXCEPTION("contract343: Left operand and result disagree about number of aux. operators.");
    }
    char trans;
    int *N, *K;
    if (transpose) {
        trans = 'T';
        N = G.params->rowtot;
        K = G.params->coltot;
    } else {
        trans = 'N';
        N = G.params->coltot;
        K = G.params->rowtot;
    }
    for (int h = 0; h < nirrep_; ++h) {
        if (b.colspi(h) != K[h]) {
            throw PSIEXCEPTION(
                "contract343: Left and right operands do not agree about the dimension of the inner index.");
        }
        if (b.colspi(h) > 0 && colspi(h) > 0) {
            global_dpd_->buf4_mat_irrep_init(&G, h);
            global_dpd_->buf4_mat_irrep_rd(&G, h);
            auto bp = b.pointer(h);
            auto rp = pointer(h);
            C_DGEMM('N', trans, b.rowspi(h), N[h], K[h], alpha, bp[0], b.colspi(h), G.matrix[h][0], G.params->coltot[h],
                    beta, rp[0], colspi(h));
        }
    }
}

DFTensor DFTensor::contract233(const Matrix& J, const DFTensor& B) {
    if (J.nirrep() != 1) {
        throw PSIEXCEPTION("contract233: Expected first argument to have no symmetry.");
    }
    auto result = DFTensor("sum_Q J(PQ) B(P|pq)", B.nQ_, B.dim2_, B.dim3_);
    auto Jp = J.pointer()[0];
    auto Jcols = J.colspi(0);
    // Sadly, we can't just make this a doublet due to symmetry...
    for (int h = 0; h < result.nirrep(); ++h) {
        if (B.colspi(h) > 0) {
            C_DGEMM('T', 'N', result.rowspi(h), result.colspi(h), B.rowspi(h), 1.0, Jp, Jcols, B.pointer(h)[0],
                    B.colspi(h), 0.0, result.pointer(h)[0], result.colspi(h));
        }
    }

    return result;
}

DFTensor DFTensor::contract123(const Matrix& Q, const Matrix& G) {
    if (Q.nirrep() != 1) {
        throw PSIEXCEPTION("contract123: Left argument must have exactly one irrep.");
    }
    if (G.symmetry()) {
        throw PSIEXCEPTION("contract123: Right argument must have trivial pont group symmetry");
    }

    auto result = DFTensor("Result", Q.colspi()[0], G.rowspi(), G.colspi());

    int offset = 0;
    for (int h = 0; h < G.nirrep(); ++h) {
        if (G.colspi(h) > 0) {
            C_DGER(Q.ncol(), G.rowspi(h) * G.colspi(h), 1.0, Q.pointer()[0], 1, G.pointer(h)[0], 1,
                   result.pointer(0)[0] + offset, result.colspi(0));
        }
        offset += G.rowspi(h) * G.colspi(h);
    }

    return result;
}

void DFTensor::add_3idx_transpose_inplace() {
    if (symmetry()) {
        throw PSIEXCEPTION("add_3idx_transpose_inplace: Tensor must be totally symmetric.");
        // In theory, ths isn't necessary, but I don't need that case.
    }
    if (dim2_ != dim3_) throw PSIEXCEPTION("DFTensor:add_3idx_transpose_inplace: Matrix must be square.");
    // Start with the totally symmetric irrep. The orbitals of each pair are of the same symmetry,
    // so treating (p, q) and (q, p) at once means we treat half the pairs within an irrep.
    int offset = 0;
    auto Mp = pointer(0);
    for (int h = 0; h < nirrep_; h++) {  // h = Irrep of first elt. in pair
        auto nh = dim2_[h];
        for (int p = 0; p < rowspi_[0]; p++) {
            for (int m = 0; m < nh; m++) {
                for (int n = 0; n <= m; n++) {
                    Mp[p][offset + m * nh + n] = Mp[p][offset + n * nh + m] =
                        Mp[p][offset + m * nh + n] + Mp[p][offset + n * nh + m];
                }
            }
        }
        offset += nh * nh;
    }
    if (colspi_[0] != offset) {
        throw PSIEXCEPTION(
            "add_3idx_transpose_inplace: Irrep 0 of Matrix isn't pairs of orbitals of appropriate symmetry from dim..");
    }
    if (nirrep_ == 1) {
        return;
    }
    // Proceed to non-totally symmetric. The orbitals of each pair are of differnet symmetires, so
    // treating (p, q) and (q, p) at once means we iterate over half the irrep pairs.
    for (int h = 1; h < nirrep_; h++) {  // h = Irrep of pair
        Mp = pointer(h);
        Dimension offsets(nirrep_);
        int offset = 0;
        for (int i = 0; i < nirrep_; i++) {
            int j = h ^ i;
            offsets[i] = offset;
            offset += dim2_[i] * dim2_[j];
        }
        for (int i = 0; i < nirrep_; i++) {  // i = Irrep of first elt. in pair
            int j = h ^ i;                      // j = Irrep of second elt. in pair
            if (j < i) {
                continue;
            }  // We already processed this pair.
            for (int p = 0; p < rowspi_[h]; p++) {
                for (int m = 0; m < dim2_[i]; m++) {
                    for (int n = 0; n < dim2_[j]; n++) {
                        Mp[p][offsets[i] + m * dim2_[j] + n] = Mp[p][offsets[j] + n * dim2_[i] + m] =
                            Mp[p][offsets[i] + m * dim2_[j] + n] + Mp[p][offsets[j] + n * dim2_[i] + m];
                    }
                }
            }
        }
        if (colspi_[h] != offset) {
            throw PSIEXCEPTION(
                "add_3idx_transpose_inplace: Matrix isn't pairs of orbitals of appropriate symmetry from dim.");
        }
    }
}

}  // namespace dct
}  // namespace psi
