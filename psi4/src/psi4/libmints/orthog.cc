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

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/orthog.h"

namespace psi {

OverlapOrthog::OverlapOrthog(OrthogMethod method, SharedMatrix overlap, double lindep_tolerance, int print)
    : orthog_method_(method), overlap_(overlap), lindep_tol_(lindep_tolerance), print_(print) {
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

void OverlapOrthog::compute_orthog_trans() {
    compute_overlap_eig();
    if (orthog_method_ == Automatic) {
        orthog_method_ = (min_S < lindep_tol_) ? Symmetric : Canonical;
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
