/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

#ifndef PSI4_GUARD_EPMP2_H
#define PSI4_GUARD_EPMP2_H

#include "psi4/libmints/wavefunction.h"

namespace psi {

// Forward declare
class DFHelper;

namespace dfep2 {

class DFEP2Wavefunction : public Wavefunction {
   protected:
    // Auxiliary basis
    std::shared_ptr<BasisSet> ribasis_;

    // Orbitals full / occ / vir / solves
    SharedMatrix AO_C_;
    SharedMatrix AO_Cocc_;
    SharedMatrix AO_Cvir_;
    SharedVector AO_eps_;
    std::vector<std::tuple<double, size_t, size_t>> orbital_order_;

    // Integrals
    std::shared_ptr<DFHelper> dfh_;

    // Options
    double conv_thresh_;
    size_t max_iter_;
    size_t debug_;
    size_t memory_doubles_;
    size_t unit_;
    size_t num_threads_;

   public:
    DFEP2Wavefunction(std::shared_ptr<Wavefunction> wfn);

    // Compute the current set of orbitals, [[h1_orb, h1_orb], [h2_orb, ...]]
    // Return [[(value1, conv1), ...], ...]
    std::vector<std::vector<std::pair<double, double>>> compute(std::vector<std::vector<size_t>> solve_orbs);
};  // End DFEP2

}  // namespace dfep2
}  // namespace psi

#endif
