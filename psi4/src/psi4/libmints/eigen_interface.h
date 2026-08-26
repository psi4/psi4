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

#pragma once

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"

#include <eigen3/Eigen/Core>
#ifdef USING_gauxc
#include <gauxc/basisset.hpp>
#endif

namespace psi {
namespace linalg {

/// Map a single-irrep Psi4 Matrix onto its underlying data.
PSI_API Eigen::Map<Eigen::MatrixXd> eigen_map(Matrix& matrix);

/// Map each irrep block of a Psi4 Matrix onto its underlying data.
PSI_API std::vector<Eigen::Map<Eigen::MatrixXd>> eigen_maps(Matrix& matrix);

Matrix matrix_from_eigen(const Eigen::MatrixXd& eigen_mat, const std::string& name = "");

// Generates the PermutationMatrix that transforms AOs from Gaussian to standard basis.
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> generate_permutation_to_cca(const BasisSet& basis);

#ifdef USING_gauxc
// converts a Psi4::BasisSet object to a GauXC::BasisSet object
template <typename T>
GauXC::BasisSet<T> to_gauxc_basisset(const BasisSet& psibasis, double basis_tol, bool force_cartesian) {
    using prim_array = typename GauXC::Shell<T>::prim_array;
    using cart_array = typename GauXC::Shell<T>::cart_array;

    auto nshell = psibasis.nshell();

    GauXC::BasisSet<T> gauxc_basisset(nshell);

    for (size_t ishell = 0; ishell != nshell; ++ishell) {
        auto psi4_shell = psibasis.shell(ishell);
        if(psi4_shell.nprimitive() > GauXC::detail::shell_nprim_max) {
            std::ostringstream oss;
            oss << "Shell has " << psi4_shell.nprimitive() <<
              " primitives but the shell datatype in GauXC is limited to " << GauXC::detail::shell_nprim_max <<
              "!\nTo work around this issue, segment your basis set!\n";
            throw std::runtime_error(oss.str());
        }

        const auto nprim = GauXC::PrimSize(psi4_shell.nprimitive());
        prim_array alpha;
        prim_array coeff;

        for (size_t iprim = 0; iprim != psi4_shell.nprimitive(); ++iprim) {
            alpha.at(iprim) = psi4_shell.exp(iprim);
            coeff.at(iprim) = psi4_shell.coef(iprim);
        }

        auto psi4_shell_center = psi4_shell.center();
        cart_array center = { psi4_shell_center[0], psi4_shell_center[1], psi4_shell_center[2] };

        gauxc_basisset[ishell] = GauXC::Shell(
            nprim,
            GauXC::AngularMomentum(psi4_shell.am()),
            (force_cartesian ? GauXC::SphericalType(false) : GauXC::SphericalType( !(psi4_shell.is_cartesian()) ) ),
            alpha,
            coeff,
            center,
            false // do not normalize shell via GauXC; it is normalized via Psi4
        );
    }

    for (auto& sh : gauxc_basisset) {
        sh.set_shell_tolerance(basis_tol);
    }

    return gauxc_basisset;
}
#endif

}  // namespace linalg
}  // namespace psi
