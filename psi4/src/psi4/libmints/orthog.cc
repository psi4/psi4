/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

OverlapOrthog::OverlapOrthog(OrthogMethod method, SharedMatrix overlap, SharedVector rsq, double lindep_tolerance,
                             double cholesky_tolerance, int print)
    : orthog_method_(method),
      overlap_(overlap),
      rsq_(rsq),
      lindep_tol_(lindep_tolerance),
      cholesky_tol_(cholesky_tolerance),
      print_(print) {
    eigval_ = nullptr;
    eigvec_ = nullptr;
}

void OverlapOrthog::compute_overlap_eig() {
    // Eigenvectors
    eigvec_ = std::make_shared<Matrix>("U", overlap_->rowspi(), overlap_->colspi());
    // Eigenvalues
    eigval_ = std::make_shared<Vector>(overlap_->colspi());
    overlap_->diagonalize(eigvec_, eigval_);

    bool min_S_initialized = false;
    for (int h = 0; h < eigval_->nirrep(); ++h) {
        for (int i = 0; i < eigval_->dim(h); i++) {
            double eval = eigval_->get(h, i);

            if (!min_S_initialized) {
                min_S = eval;
                min_S_initialized = true;
            } else if (min_S > eval) {
                min_S = eval;
            }
        }
    }
    outfile->Printf("  Minimum eigenvalue in the overlap matrix is %14.10E.\n", min_S);
}

void OverlapOrthog::compute_inverse() {
    Xinv_ = std::make_shared<Matrix>("Orthogonal Inverse Transformation", X_->rowspi(), X_->colspi());
    Xinv_->gemm(false, false, 1.0, overlap_, X_, 0.0);
}

SharedMatrix OverlapOrthog::overlap_inverse() {
    auto Sinv = std::make_shared<Matrix>("Inverse Transformation", X_->rowspi(), X_->rowspi());
    Sinv->gemm(false, true, 1.0, X_, X_, 0.0);
    return Sinv;
}

void OverlapOrthog::compute_symmetric_orthog() {
    if (!eigval_) compute_overlap_eig();
    if (min_S < lindep_tol_) {
        outfile->Printf("WARNING: smallest overlap eigenvalue %e is smaller than S_TOLERANCE!\n", min_S);
    }
    const Dimension& nbf = eigval_->dimpi();

    SharedVector eval(std::make_shared<Vector>(nbf));
    eval->copy(*eigval_);
    for (int h = 0; h < eval->nirrep(); ++h) {
        for (int i = 0; i < eval->dim(h); i++) {
            double val = 1.0 / sqrt(eval->get(h, i));
            eval->set(h, i, val);
        }
    }

    auto eigtemp2 = std::make_shared<Matrix>("Scratch matrix 1", nbf, nbf);
    eigtemp2->set_diagonal(eval);
    auto eigtemp = std::make_shared<Matrix>("Scratch matrix 2", nbf, nbf);
    eigtemp->gemm(false, true, 1.0, eigtemp2, eigvec_, 0.0);

    X_ = std::make_shared<Matrix>("X (Symmetric Orthogonalization)", nbf, nbf);
    X_->gemm(false, false, 1.0, eigvec_, eigtemp, 0.0);
}

void OverlapOrthog::compute_canonical_orthog() {
    if (!eigval_) compute_overlap_eig();
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

std::vector<std::vector<int>> OverlapOrthog::sort_indices() const {
    std::vector<std::vector<int>> order(rsq_->nirrep());
    for (int h = 0; h < rsq_->nirrep(); h++) {
        // initialize ordering as 0, 1, 2, ..., n-1
        order[h].resize(rsq_->dim(h));
        iota(order[h].begin(), order[h].end(), 0);

        // and then sort the indices by rsq
        std::sort(order[h].begin(), order[h].end(),
                  [this, h](int i1, int i2) { return rsq_->get(h, i1) < rsq_->get(h, i2); });
    }

    return order;
}

void OverlapOrthog::compute_partial_cholesky_orthog() {
    // Original dimensions
    const Dimension& nbf = overlap_->rowspi();

    // Order basis functions from tight to diffuse
    std::vector<std::vector<int>> order = sort_indices();
    // Reorder overlap matrix
    auto Stmp = std::make_shared<Matrix>(nbf, nbf);
    for (size_t h = 0; h < order.size(); h++) {
        for (size_t i = 0; i < order[h].size(); i++) {
            for (size_t j = 0; j <= i; j++) {
                Stmp->set(h, i, j, overlap_->get(h, order[h][i], order[h][j]));
                Stmp->set(h, j, i, overlap_->get(h, order[h][i], order[h][j]));
            }
        }
    }
    // Do pivoted Cholesky to find the important basis functions
    std::vector<std::vector<int>> pivots;
    Stmp->partial_cholesky_factorize_pivot(cholesky_tol_, true, pivots);
    // Rewrite the pivot indices in terms of the original indexing
    for (size_t h = 0; h < pivots.size(); h++)
        for (size_t i = 0; i < pivots[h].size(); i++) pivots[h][i] = order[h][pivots[h][i]];

    // Size of Cholesky basis
    Dimension nchol(nbf);
    for (int h = 0; h < overlap_->nirrep(); h++) {
        nchol[h] = pivots[h].size();
        outfile->Printf("  Cholesky: irrep %d, %d of %d possible MOs eliminated.\n", h, nbf[h] - nchol[h], nbf[h]);
    }
    if (nchol.sum() != nbf.sum() && print_)
        outfile->Printf("  Cholesky: overall, %d of %d possible MOs eliminated.\n\n", nbf.sum() - nchol.sum(),
                        nbf.sum());

    // Copy over data
    auto Ssub = std::make_shared<Matrix>(nchol, nchol);
    for (int h = 0; h < overlap_->nirrep(); h++)
        for (int m = 0; m < (int)pivots[h].size(); m++)
            for (int n = 0; n < (int)pivots[h].size(); n++)
                Ssub->set(h, m, n, overlap_->get(h, pivots[h][m], pivots[h][n]));

    // Switch overlap matrix to use
    auto overlap0 = overlap_;
    overlap_ = Ssub;
    // Reset eigendecomposition
    eigvec_ = nullptr;
    eigval_ = nullptr;

    // Compute orthogonalization as usual in submatrix
    if (print_) outfile->Printf("    Proceeding with conventional orthogonalization in reduced basis.\n");
    orthog_method_ = Automatic;
    compute_orthog_trans();

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
    overlap_ = overlap0;
    // Get orthogonal transformation
    X_ = padX;
}

void OverlapOrthog::compute_orthog_trans() {
    compute_overlap_eig();
    if (orthog_method_ == Automatic) {
        orthog_method_ = (min_S < lindep_tol_) ? Canonical : Symmetric;
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
            if (print_) outfile->Printf("    Using partial Cholesky orthogonalization.\n");
            compute_partial_cholesky_orthog();
            break;
        default:
            throw PSIEXCEPTION("OverlapOrthog::compute_orthog_trans: bad value.");
    }

    // Compute inverse transformation, too.
    compute_inverse();
}

SharedMatrix OverlapOrthog::basis_to_orthog_basis() {
    if (!X_) compute_orthog_trans();

    return X_;
}

SharedMatrix OverlapOrthog::basis_to_orthog_basis_inverse() {
    if (!X_) compute_orthog_trans();
    if (!Xinv_) compute_inverse();

    return Xinv_;
}

Dimension OverlapOrthog::dim() { return X_->rowdim(); }

Dimension OverlapOrthog::orthog_dim() {
    if (!X_) compute_orthog_trans();

    return X_->coldim();
}

int OverlapOrthog::nlindep() {
    if (!X_) compute_orthog_trans();

    return X_->colspi().sum();
}

int OverlapOrthog::nlindep(int h) {
    if (!X_) compute_orthog_trans();

    return X_->coldim(h);
}

}  // namespace psi
