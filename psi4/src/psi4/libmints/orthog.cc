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

#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/orthog.h"

namespace psi {

BasisSetOrthogonalization::BasisSetOrthogonalization(OrthogonalizationMethod method, SharedMatrix overlap,
                                                     double lindep_tolerance, double cholesky_tolerance, int print)
    : orthog_method_(method),
      overlap_(overlap),
      lindep_tol_(lindep_tolerance),
      cholesky_tol_(cholesky_tolerance),
      print_(print) {
    eigval_ = nullptr;
    eigvec_ = nullptr;
    normalization_ = nullptr;
    normalized_overlap_ = nullptr;
}

void BasisSetOrthogonalization::normalize() {
    // Compute normalization coefficients
    if (normalization_ != nullptr)
        throw PSIEXCEPTION("BasisSetOrthogonalization::normalize: normalization_ should be nullptr");
    normalization_ = std::make_shared<Vector>(overlap_->rowspi());
    normalization_->set_name("SO normalization factors");
    for (int h = 0; h < overlap_->nirrep(); h++)
        for (int i = 0; i < overlap_->rowdim(h); i++) normalization_->set(h, i, 1.0 / sqrt(overlap_->get(h, i, i)));
    if (print_ > 3) normalization_->print();

    // Normalize overlap matrix
    if (normalized_overlap_ != nullptr)
        throw PSIEXCEPTION("BasisSetOrthogonalization::normalize: normalized_overlap_ should be nullptr");
    normalized_overlap_ = std::make_shared<Matrix>(overlap_->rowspi(), overlap_->colspi());
    for (int h = 0; h < overlap_->nirrep(); h++)
        for (int i = 0; i < overlap_->rowdim(h); i++)
            for (int j = 0; j < overlap_->coldim(h); j++)
                normalized_overlap_->set(
                    h, i, j, overlap_->get(h, i, j) * normalization_->get(h, i) * normalization_->get(h, j));
}

void BasisSetOrthogonalization::unroll_normalization() {
    if (!normalization_)
        throw PSIEXCEPTION("BasisSetOrthogonalization::unroll_normalization: normalization has not been yet computed.");
    if (!X_) throw PSIEXCEPTION("BasisSetOrthogonalization::unroll_normalization: X has not been yet computed.");

    // Scaling the rows by the normalization coefficients leads to
    // X^T*S*X = 1, since S is thereby transformed to the normalized
    // basis used to form X.
    for (int h = 0; h < X_->nirrep(); h++)
        for (int i = 0; i < X_->rowdim(h); i++) X_->scale_row(h, i, normalization_->get(h, i));
}

void BasisSetOrthogonalization::compute_overlap_eig() {
    if (!normalized_overlap_)
        throw PSIEXCEPTION(
            "BasisSetOrthogonalization::compute_overlap_eig: normalized overlap has not yet been computed.");

    // Eigenvectors
    eigvec_ = std::make_shared<Matrix>("U", normalized_overlap_->rowspi(), normalized_overlap_->colspi());
    // Eigenvalues
    eigval_ = std::make_shared<Vector>(normalized_overlap_->colspi());
    normalized_overlap_->diagonalize(eigvec_, eigval_);

    // Find minimum eigenvalue
    bool min_S_initialized = false;
    for (int h = 0; h < eigval_->nirrep(); ++h) {
        for (int i = 0; i < eigval_->dim(h); i++) {
            double eval = eigval_->get(h, i);
            if (!min_S_initialized) {
                min_S_ = eval;
                min_S_initialized = true;
            } else if (min_S_ > eval) {
                min_S_ = eval;
            }
        }
    }
    if (print_) outfile->Printf("  Minimum eigenvalue in the overlap matrix is %14.10E.\n", min_S_);

    // Find reciprocal condition number
    rcond_ = DBL_MAX;
    for (int h = 0; h < eigval_->nirrep(); ++h) {
        if (!eigval_->dim(h)) continue;
        // Min and max eigenvalue in symmetry block
        double e_min, e_max;
        e_min = e_max = eigval_->get(h, 0);
        for (int i = 0; i < eigval_->dim(h); i++) {
            double eval = eigval_->get(h, i);
            e_min = std::min(e_min, eval);
            e_max = std::max(e_max, eval);
        }
        // Reciprocal condition number
        rcond_ = std::min(rcond_, e_min / e_max);
    }
    if (print_) outfile->Printf("  Reciprocal condition number of the overlap matrix is %14.10E.\n", rcond_);
}

void BasisSetOrthogonalization::compute_inverse() {
    Xinv_ = std::make_shared<Matrix>("Orthogonal Inverse Transformation", X_->rowspi(), X_->colspi());
    Xinv_->gemm(false, false, 1.0, overlap_, X_, 0.0);
}

SharedMatrix BasisSetOrthogonalization::overlap_inverse() {
    auto Sinv = std::make_shared<Matrix>("Inverse overlap matrix", X_->rowspi(), X_->rowspi());
    Sinv->gemm(false, true, 1.0, X_, X_, 0.0);
    return Sinv;
}

void BasisSetOrthogonalization::compute_symmetric_orthog() {
    /*
      Susi Lehtola Jan 16 2020

      SO basis functions aren't normalized, which causes X to be
      asymmetric at the end (see unroll_normalization function). It
      might be possible to achieve a symmetric X even with the use of
      balancing i.e. renormalization of the basis set, but this would
      require tracking the dependence of the eigenvectors and
      eigenvalues on the scaling transform.

      If one omits the normalization with the commented code below, X
      becomes symmetric. However, the symmetricity of X or the lack of
      it makes no difference at the end, since the orbitals anyway arise
      from diagonalizing the Fock matrix or orbital optimization, and
      any code in Psi4 should also work with canonical orthogonalization
      which is asymmetric by force.
    */

    /*
      // Bypass re-normalization of the basis set
      normalized_overlap_ = overlap_;
      for (int h = 0; h < overlap_->nirrep(); h++)
        for (int i = 0; i < overlap_->rowdim(h); i++) normalization_->set(h, i, 1.0);
      compute_overlap_eig();
    */

    if (!eigval_) compute_overlap_eig();
    if (min_S_ < lindep_tol_) {
        outfile->Printf("WARNING: smallest overlap eigenvalue %e is smaller than S_TOLERANCE!\n", min_S_);
    }
    const Dimension& nbf = eigval_->dimpi();
    int nirrep = eigval_->nirrep();

    // Symmetric orthogonalization obtained as U s^{-1/2} U^t
    auto Us = std::make_shared<Matrix>("Half-transformed matrix Us", nbf, nbf);
    Us->copy(eigvec_);
    for (int h = 0; h < nirrep; ++h) {
        for (int i = 0; i < nbf[h]; ++i) {
            Us->scale_column(h, i, 1.0 / sqrt(eigval_->get(h, i)));
        }
    }
    X_ = std::make_shared<Matrix>("X (Symmetric Orthogonalization)", nbf, nbf);
    // X = U s U^t
    X_->gemm(false, true, 1.0, Us, eigvec_, 0.0);
}

void BasisSetOrthogonalization::compute_canonical_orthog() {
    if (!eigval_) compute_overlap_eig();
    if (rcond_ <= DBL_EPSILON)
        outfile->Printf(
            "WARNING: overlap condition number is close to machine precision %14.10E, indicating a pathologically "
            "overcomplete basis. Think about switching to partial Cholesky.\n",
            DBL_EPSILON);

    const Dimension& nbf = eigval_->dimpi();
    int nirrep = eigval_->nirrep();

    // Count number of orbitals
    Dimension nmo(nbf);
    for (int h = 0; h < nirrep; ++h) {
        nmo[h] = 0;
        for (int i = 0; i < nbf[h]; ++i) {
            if (eigval_->get(h, i) > lindep_tol_) {
                nmo[h]++;
            }
        }
        if (print_ > 2) outfile->Printf("  Irrep %d, %d of %d possible MOs eliminated.\n", h, nbf[h] - nmo[h], nbf[h]);
    }
    if (nmo.sum() != nbf.sum() && print_)
        outfile->Printf("  Overall, %d of %d possible MOs eliminated.\n\n", nbf.sum() - nmo.sum(), nbf.sum());

    // Form vectors
    X_ = std::make_shared<Matrix>("X (Canonical Orthogonalization)", nbf, nmo);
    for (int h = 0; h < nirrep; ++h) {
        for (int mo = 0; mo < nmo[h]; ++mo) {
            // Columns are in increasing overlap eigenvalue
            int col = nbf[h] - nmo[h] + mo;
            double norm = 1.0 / sqrt(eigval_->get(h, col));
            for (int j = 0; j < nbf[h]; j++) X_->set(h, j, mo, norm * eigvec_->get(h, j, col));
        }
    }
}

std::vector<std::vector<int>> BasisSetOrthogonalization::sort_indices() const {
    std::vector<std::vector<int>> order(normalized_overlap_->nirrep());
    for (int h = 0; h < normalized_overlap_->nirrep(); h++) {
        // initialize ordering as 0, 1, 2, ..., n-1
        order[h].resize(normalized_overlap_->coldim(h));
        iota(order[h].begin(), order[h].end(), 0);

        // Collect off-diagonal overlap
        std::vector<double> od(normalized_overlap_->coldim(h));
        for (size_t i = 0; i < od.size(); i++) {
            double sum = 0.0;
            for (size_t j = 0; j < od.size(); j++) {
                if (i == j) continue;
                sum += std::abs(normalized_overlap_->get(h, i, j));
            }
            od[i] = sum;
        }

        // and then sort the indices by increasing off-diagonal overlap
        std::stable_sort(order[h].begin(), order[h].end(), [&od](int i1, int i2) { return od[i1] < od[i2]; });
    }

    return order;
}

void BasisSetOrthogonalization::compute_partial_cholesky_orthog() {
    // Original dimensions
    const Dimension& nbf = normalized_overlap_->rowspi();

    // Get the best initial order for the basis functions
    std::vector<std::vector<int>> order = sort_indices();
    // Reorder overlap matrix
    auto Stmp = std::make_shared<Matrix>(nbf, nbf);
    for (size_t h = 0; h < order.size(); h++) {
        for (size_t i = 0; i < order[h].size(); i++) {
            for (size_t j = 0; j <= i; j++) {
                Stmp->set(h, i, j, normalized_overlap_->get(h, order[h][i], order[h][j]));
                Stmp->set(h, j, i, normalized_overlap_->get(h, order[h][i], order[h][j]));
            }
        }
    }

    // Do pivoted Cholesky to find the important basis functions
    std::vector<std::vector<int>> pivots;
    Stmp->pivoted_cholesky(cholesky_tol_, pivots);
    // and rewrite the pivot indices in terms of the original indexing
    for (size_t h = 0; h < pivots.size(); h++)
        for (size_t i = 0; i < pivots[h].size(); i++) pivots[h][i] = order[h][pivots[h][i]];
    if (print_ > 2) {
        outfile->Printf("    Cholesky pivot functions:\n");
        for (size_t h = 0; h < pivots.size(); h++)
            for (size_t i = 0; i < pivots[h].size(); i++)
                outfile->Printf("    Symmetry %u function %u: SO basis function %i\n", (unsigned int)h + 1,
                                (unsigned int)i + 1, pivots[h][i] + 1);
    }

    // Size of Cholesky basis
    Dimension nchol(nbf);
    for (int h = 0; h < normalized_overlap_->nirrep(); h++) {
        nchol[h] = pivots[h].size();
        if (print_ > 2)
            outfile->Printf("  Cholesky: irrep %d, %d of %d possible AOs eliminated.\n", h, nbf[h] - nchol[h], nbf[h]);
    }
    if (nchol.sum() != nbf.sum() && print_)
        outfile->Printf("  Cholesky: overall, %d of %d possible AOs eliminated.\n\n", nbf.sum() - nchol.sum(),
                        nbf.sum());

    // Copy over data
    auto Ssub = std::make_shared<Matrix>(nchol, nchol);
    for (int h = 0; h < normalized_overlap_->nirrep(); h++)
        for (int m = 0; m < (int)pivots[h].size(); m++)
            for (int n = 0; n < (int)pivots[h].size(); n++)
                Ssub->set(h, m, n, normalized_overlap_->get(h, pivots[h][m], pivots[h][n]));

    // Switch overlap matrix to use
    auto normalized_overlap0 = normalized_overlap_;
    normalized_overlap_ = Ssub;
    // Reset eigendecomposition
    eigvec_ = nullptr;
    eigval_ = nullptr;

    // Compute orthogonalization as usual in submatrix
    if (print_) outfile->Printf("    Proceeding with canonical orthogonalization in reduced basis.\n");
    compute_canonical_orthog();

    // Find out number of vectors in each symmetry block
    const Dimension& nmo(X_->colspi());

    // Padded matrices
    SharedMatrix padX = std::make_shared<Matrix>("Orthogonal Transformation", nbf, nmo);
    padX->zero();
    for (int h = 0; h < X_->nirrep(); h++)
        for (int n = 0; n < nmo[h]; n++)
            for (int m = 0; m < (int)pivots[h].size(); m++) padX->set(h, pivots[h][m], n, X_->get(h, m, n));

    orthog_method_ = PartialCholesky;
    // Restore original overlap matrix
    normalized_overlap_ = normalized_overlap0;
    // Get orthogonal transformation
    X_ = padX;
}

void BasisSetOrthogonalization::compute_orthog_trans() {
    // Normalize basis
    normalize();

    // Determine what method to use
    if (orthog_method_ == Automatic) {
        compute_overlap_eig();
        if (rcond_ <= DBL_EPSILON || !std::isnormal(rcond_))
            orthog_method_ = PartialCholesky;
        else if (min_S_ < lindep_tol_)
            orthog_method_ = Canonical;
        else
            orthog_method_ = Symmetric;
    }

    switch (orthog_method_) {
        case Symmetric:
            if (print_) outfile->Printf("    Using symmetric orthogonalization.\n");
            compute_symmetric_orthog();
            break;
        case Canonical:
            if (print_) outfile->Printf("    Using canonical orthogonalization.\n");
            compute_canonical_orthog();
            break;
        case PartialCholesky:
            if (print_) outfile->Printf("    Using partial Cholesky orthogonalization (doi:10.1063/1.5139948, doi:10.1103/PhysRevA.101.032504).\n");
            compute_partial_cholesky_orthog();
            break;
        default:
            throw PSIEXCEPTION("BasisSetOrthogonalization::compute_orthog_trans: bad value.");
    }

    // Include basis function normalization in X
    unroll_normalization();
    // Compute the inverse transformation, too.
    compute_inverse();
    // Check that the computed basis truly is orthonormal
    check_orth();
}

void BasisSetOrthogonalization::check_orth() {
    // Check that orbitals are orthonormal
    const Dimension& nbf(X_->rowspi());
    const Dimension& nmo(X_->colspi());
    auto SX = std::make_shared<Matrix>("SX", nbf, nmo);
    SX->gemm(false, false, 1.0, overlap_, X_, 0.0);
    auto XtSX = std::make_shared<Matrix>("MO overlap", nmo, nmo);
    XtSX->gemm(true, false, 1.0, X_, SX, 0.0);

    if (print_ > 3) XtSX->print();

    // Remove unity from diagonal
    for (int h = 0; h < X_->nirrep(); h++)
        for (int n = 0; n < nmo[h]; n++) XtSX->set(h, n, n, XtSX->get(h, n, n) - 1.0);

    // Compute norm
    double norm = 0.0;
    for (int h = 0; h < X_->nirrep(); h++)
        for (int m = 0; m < nmo[h]; m++)
            for (int n = 0; n < nmo[h]; n++) norm += std::pow(XtSX->get(h, m, n), 2);

    if (print_ > 2) outfile->Printf("  MO non-orthonormality %e\n", norm);

    if (norm >= 1e-10) throw PSIEXCEPTION("BasisSetOrthogonalization::check_orth: orbitals are not orthonormal");
}

SharedMatrix BasisSetOrthogonalization::basis_to_orthog_basis() {
    if (!X_) compute_orthog_trans();

    return X_;
}

SharedMatrix BasisSetOrthogonalization::basis_to_orthog_basis_inverse() {
    if (!X_) compute_orthog_trans();
    if (!Xinv_) compute_inverse();

    return Xinv_;
}

Dimension BasisSetOrthogonalization::dim() { return X_->rowdim(); }

Dimension BasisSetOrthogonalization::orthog_dim() {
    if (!X_) compute_orthog_trans();

    return X_->coldim();
}

int BasisSetOrthogonalization::nlindep() {
    if (!X_) compute_orthog_trans();

    return X_->colspi().sum();
}

int BasisSetOrthogonalization::nlindep(int h) {
    if (!X_) compute_orthog_trans();

    return X_->coldim(h);
}

}  // namespace psi
