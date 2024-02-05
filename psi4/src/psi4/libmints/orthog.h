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

#ifndef ORTHOG_H
#define ORTHOG_H

#include "typedefs.h"
#include "psi4/libmints/dimension.h"

namespace psi {

/*! \ingroup MINTS
 *  \class BasisSetOrthogonalization
 *  \brief Implements methods for orthogonalizing basis sets.
 */
class PSI_API BasisSetOrthogonalization {
   public:
    /// An enum for the types of orthogonalization.
    enum OrthogonalizationMethod { Symmetric, Canonical, PartialCholesky, Automatic };

   private:
    int print_;

    /// Matrix to decompose
    SharedMatrix overlap_;
    /// Normalized version of the input matrix
    SharedMatrix normalized_overlap_;
    /// Normalization coefficients
    SharedVector normalization_;
    /// its eigenvectors
    SharedMatrix eigvec_;
    /// and eigenvalues
    SharedVector eigval_;

    /** The tolerance for linearly independent basis functions.
     * The interpretation depends on the orthogonalization
     * method.
     */
    double lindep_tol_;
    /// The Cholesky decomposition threshold
    double cholesky_tol_;
    /// The orthogonalization method
    OrthogonalizationMethod orthog_method_;
    /// The orthogonalizing matrix
    SharedMatrix X_;
    /// ... and its inverse
    SharedMatrix Xinv_;

    /// Smallest eigenvalue
    double min_S_;
    /// Reciprocal condition number
    double rcond_;
    /// @}

    /** Normalizes the basis set. This is important for the partial
        Cholesky decomposition, as otherwise the functions with most
        overlap i.e. the most diffuse ones get handled first by the
        algorithm, which is exactly the wrong way around.
    */
    void normalize();
    /// Once X has been formed, unroll the normalization into X
    void unroll_normalization();

    /// Given X, compute Xinv = S*X
    void compute_inverse();
    /// Compute eigendecomposition
    void compute_overlap_eig();
    /// Symmetric orthogonalization
    void compute_symmetric_orthog();
    /// Canonical orthogonalization
    void compute_canonical_orthog();
    /// Partial Cholesky orthogonalization
    void compute_partial_cholesky_orthog();
    /// Driver routine
    void compute_orthog_trans();
    /// Check computed basis is orthonormal
    void check_orth();
    /// Sort the basis functions from tight to diffuse
    std::vector<std::vector<int>> sort_indices() const;

   public:
    BasisSetOrthogonalization(OrthogonalizationMethod method, SharedMatrix overlap, double lindep_tolerance,
                              double cholesky_tolerance, int print = 0);

    OrthogonalizationMethod orthog_method() const { return orthog_method_; }

    double lindep_tol() const { return lindep_tol_; }

    /** Returns a matrix which does the requested transform
        from a basis to an orthogonal basis.  This could be
        either the symmetric or canonical orthogonalization
        matrix.  The row dimension is the basisdimension and
        the column dimension is orthogonal basis dimension.An
        operator \f$O\f$ in the orthogonal basis is given by
        \f$X OX^T\f$ where \f$X\f$ is the return value of this
        function.
    */
    SharedMatrix basis_to_orthog_basis();

    /** Returns the inverse of the transformation returned by
     * basis_to_orthog_basis().
     */
    SharedMatrix basis_to_orthog_basis_inverse();

    /// Return an $S^{-1}$.
    SharedMatrix overlap_inverse();

    /// Number of basis functions
    Dimension dim();
    /// Number of orthogonal functions
    Dimension orthog_dim();

    /// Number of independent functions
    int nlindep();
    /// Number of independent functions in symmetry block h
    int nlindep(int h);
};

}  // namespace psi

#endif  // ORTHOG_H
