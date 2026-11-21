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

#ifndef _psi_src_lib_libmints_batchedpotential_h_
#define _psi_src_lib_libmints_batchedpotential_h_

#include <array>
#include <memory>
#include <vector>

#include "psi4/pragma.h"

namespace libint2 {
template <typename Real>
class FmEval_Chebyshev7;
class Shell;
}  // namespace libint2

namespace psi {

class BasisSet;
class Matrix;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class BatchedPotentialInt
 *  \brief Potential integrals (mu| 1/|r - C_k| |nu) of one shell pair over a batch
 *         of points C_k, with the point index innermost.
 *
 *  These are the integrals seminumerical exchange (COSX), continuum solvation, and
 *  electrostatic potentials on grids need at many points. Evaluating a shell pair
 *  over all the points at once, rather than one point per call, lets the work that
 *  depends only on the shell pair be done once and the rest vectorize over the
 *  points.
 *
 *  The integrals are computed with the McMurchie-Davidson scheme: per primitive pair,
 *  the Hermite expansion coefficients (independent of the points) are formed once,
 *  and the Hermite Coulomb integrals R_tuv are formed for a chunk of points at a time.
 *  They are positive for positive overlap densities: the potential of a unit positive
 *  charge acting on a positive charge distribution, as ElectrostaticInt computes them.
 */
class PSI_API BatchedPotentialInt {
   public:
    BatchedPotentialInt(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2);
    ~BatchedPotentialInt();

    /// Compute the integrals of bra shell s1 (of bs1) and ket shell s2 (of bs2) at the
    /// npoints points (x[k], y[k], z[k]) [Bohr]. The result is in buffer(), as
    /// [i1][i2][k]: function i1 of s1, function i2 of s2 (in psi4's order, spherical or
    /// cartesian as each basis set asks), and point k, which runs fastest.
    void compute_shell(int s1, int s2, size_t npoints, const double* x, const double* y, const double* z);

    /// The integrals the last compute_shell computed.
    const double* buffer() const { return out_.data(); }

    std::shared_ptr<BasisSet> basis1() const { return bs1_; }
    std::shared_ptr<BasisSet> basis2() const { return bs2_; }

   private:
    /// F_n(T[k]) into F[n * nk + k], n <= L, for the nk arguments T; inv_T is nk doubles of work space
    void boys(int L, size_t nk, const double* T, double* F, double* inv_T);
    /// Cartesian integrals into cart_, by the generated Obara-Saika VRR or HRR kernels (angular momenta up to
    /// snk::max_am)
    void compute_cart_vrr(const libint2::Shell& sh1, const libint2::Shell& sh2, size_t npoints, const double* x,
                          const double* y, const double* z);
    /// Cartesian integrals into cart_, by the McMurchie-Davidson scheme (any angular momenta)
    void compute_cart_md(const libint2::Shell& sh1, const libint2::Shell& sh2, size_t npoints, const double* x,
                         const double* y, const double* z);

    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    /// Highest total angular momentum of a shell pair
    int max_L_;
    std::shared_ptr<const libint2::FmEval_Chebyshev7<double>> fm_eval_;
    /// Cartesian -> spherical transforms by angular momentum
    std::vector<SphericalTransform> sph_trans_;
    /// Cartesian components (lx, ly, lz) by angular momentum, in psi4's (CCA) order
    std::vector<std::vector<std::array<int, 3>>> comps_;

    /// Hermite expansion coefficients of a primitive pair, per direction
    std::vector<double> Ex_, Ey_, Ez_;
    /// Boys function values of a chunk of points, [n][k], then T_k and 1/T_k
    std::vector<double> boys_;
    /// Hermite Coulomb integrals of a chunk of points, [n][t][u][v][k]
    std::vector<double> R_;
    /// Cartesian integrals of the shell pair, [c1][c2][k], and a half-transformed copy
    std::vector<double> cart_;
    std::vector<double> half_;
    /// The result, [i1][i2][k]
    std::vector<double> out_;

    /// Work space of compute_cart_vrr: per significant primitive pair p, P, P - A, and (2 pi / p) U_p; P - A by
    /// direction; and for a chunk of points, P - G, T, the seeds, and the kernel's output
    std::vector<double> prim_p_, prim_P_, prim_PA_, prim_U_;
    std::vector<double> vrr_PA_, vrr_PG_, vrr_T_, vrr_seed_, vrr_out_, hrr_base_;
    /// Work space of boys: the interpolated Boys function values, [n][k], then exp(-T)
    std::vector<double> table_F_;
};

/// The electrostatic potential integrals contracted with a density matrix at many points:
/// V[k] += sum_{mu nu} D_{mu nu} (mu| 1/|r - C_k| |nu), for mu in bs1 and nu in bs2.
/// D is nbf(bs1) x nbf(bs2), points is npoints x 3 [Bohr], and V has npoints elements.
/// Threaded over shell pairs.
PSI_API void batched_potential_contract_density(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                                const Matrix& D, const Matrix& points, double* V);

/// The potential matrix of point charges q[k] at many points:
/// F_{mu nu} += sum_k q[k] (mu| 1/|r - C_k| |nu), for mu in bs1 and nu in bs2.
/// F is nbf(bs1) x nbf(bs2), and points is npoints x 3 [Bohr]. Threaded over shell pairs.
PSI_API void batched_potential_contract_charges(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                                const double* q, const Matrix& points, Matrix& F);

}  // namespace psi

#endif
