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

#ifndef _psi_src_lib_libmints_sox2c1e_h_
#define _psi_src_lib_libmints_sox2c1e_h_

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/dimension.h"

#include <Einsums/Config.hpp>
#include <Einsums/Tensor.hpp>

using ComplexMatrix = einsums::BlockTensor<std::complex<double>, 2>;
using SharedComplexMatrix = std::shared_ptr<ComplexMatrix>;

namespace psi {

class IntegralFactory;
class MatrixFactory;
class BasisSet;

class PSI_API SOX2C1e {
   public:
    SOX2C1e(std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    ~SOX2C1e();

    void compute(SharedMatrix S, SharedMatrix T, SharedMatrix V, SharedMatrix Hx, SharedMatrix Hy, SharedMatrix Hz);

   private:
    /// Compute the S, T, V, and W integrals
    void compute_integrals(std::shared_ptr<IntegralFactory>, std::shared_ptr<MatrixFactory>);
    void form_dirac_hamiltonian(SharedComplexMatrix, SharedComplexMatrix);
    void form_orth(SharedComplexMatrix);
    void form_X(ComplexMatrix const&, SharedComplexMatrix);
    void form_Stilde(ComplexMatrix const&, SharedComplexMatrix);
    void form_R(ComplexMatrix const&, ComplexMatrix const&, SharedComplexMatrix);

    /// The contracted AO basis
    std::shared_ptr<BasisSet> contractedBasis_;
    /// The uncontracted AO basis used to solve the Dirac equation
    std::shared_ptr<BasisSet> uncBasis_;
    /// Dimension of the orbital basis and 2/4-component bases
    Dimension nsopi_;
    Dimension nsopi2c_;
    Dimension nsopi4c_;

    // 1e Integrals in 2-component complex form
    SharedComplexMatrix overlap_;
    SharedComplexMatrix kinetic_;
    SharedComplexMatrix nuclear_;
    // relativistic potential W (\sigma\cdot p)V(\sigma\cdot p)
    SharedComplexMatrix rel_pot_;

};

}  // namespace psi

#endif
