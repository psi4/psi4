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

#ifndef psi4_libscf_solver_cghf_h
#define psi4_libscf_solver_cghf_h

#include "psi4/libmints/complexwavefunction.h"
#include "psi4/libpsio/psio.hpp"

#include <vector>

namespace psi {
class BasisSet;
class SuperFunctional;
namespace scf {

class CGHF : public ComplexWavefunction {
   protected:
    /// The DFT functional definition
    std::shared_ptr<SuperFunctional> functional_;

    /// Basis list for SAD
    std::vector<std::shared_ptr<BasisSet>> sad_basissets_;
    std::vector<std::shared_ptr<BasisSet>> sad_fitting_basissets_;

    void common_init();

   public:
    CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> functional);
    CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> functional, Options& options,
         std::shared_ptr<PSIO> psio);
    ~CGHF() override;

    /// SAD information
    void set_sad_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) { sad_basissets_ = basis_vec; }
    void set_sad_fitting_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) {
        sad_fitting_basissets_ = basis_vec;
    }
};

}  // namespace scf
}  // namespace psi

#endif
