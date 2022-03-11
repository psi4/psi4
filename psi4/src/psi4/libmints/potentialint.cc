/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "psi4/libmints/potentialint.h"
#include "psi4/libpsi4util/process.h"
#include "libint2/engine.h"

namespace psi {

PCMPotentialInt::PCMPotentialInt(std::vector<SphericalTransform>& trans, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> /* bs2 */, int /* deriv */)
    : PotentialInt(trans, bs1, bs1) {

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    int nthreads = Process::environment.get_n_threads();
    for(int thread=0; thread < nthreads; ++thread) {
        engines_.push_back(std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 0));
    }
}

}  // namespace psi
