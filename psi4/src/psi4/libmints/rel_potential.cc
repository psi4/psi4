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

#include "psi4/libmints/rel_potential.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"

#include <libint2/shell.h>
#include <libint2/engine.h>

#include <algorithm>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define RELVDEBUG 0

;
using namespace psi;

// Initialize potential_recur_ to +1 basis set angular momentum
RelPotentialInt::RelPotentialInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2, int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv) {
    if (deriv > 0) {
        throw PSIEXCEPTION("RelPotentialInt: deriv > 0 is not supported.");
    }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());
    engine2_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 2);

    buffer_ = new double[INT_NCART(max_am) * INT_NCART(max_am)];
    buffers_.resize(1);
    buffers_[0] = buffer_;
}

RelPotentialInt::~RelPotentialInt() {
    delete[] buffer_;
}

/*
 * This code was originally written by Prakash in the Evangelista lab, but modified
 * by Andy Simmonett to use Libint2.  We need to compute  <p mu |1/r-C | p nu > where
 * p is the del operator; these are easily obtained from second derivative integrals.
*/
void RelPotentialInt::compute_pair(const libint2::Shell& s1, const libint2::Shell& s2) {
    size_t size = s1.size() * s2.size();
    ::memset(buffer_, 0, size * sizeof(double));

    // If we only add one center, C, at a time, Libint2 is more memory efficient because
    // it creates buffers corresponding to each external point charge..  It also yields
    // predictable ordering of integral buffers, which are organized as follows...
    //
    //   0    1    2    3    4    5    6    7    8    9   10
    // AxAx AxAy AxAz AxBx AxBy AxBz AxCx AxCy AxCz AyAy AyAz
    //  11   12   13   14   15   16   17   18   19   20  ....
    // AyBx AyBy AyBz AyCx AyCy AyCz AzAz AzBx AzBy AzBz ....
    //
    const auto &results = engine2_->results();
    for (int A = 0; A < bs1_->molecule()->natom(); A++) {
        const auto &mol = *bs1_->molecule();
        // Setup the initial field of partial charges
        engine2_->set_params(std::vector<std::pair<double, std::array<double, 3>>>{{mol.Z(A),{mol.x(A), mol.y(A), mol.z(A)}}});
        engine2_->compute(s1, s2);
        // Add AxBx
        std::transform(buffer_, buffer_+size, results[3], buffer_, std::plus<>{});
        // Add AyBy
        std::transform(buffer_, buffer_+size, results[12], buffer_, std::plus<>{});
        // Add AzBz
        std::transform(buffer_, buffer_+size, results[20], buffer_, std::plus<>{});
    }
}

RelPotentialSOInt::RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>& aoint,
                                     const std::shared_ptr<IntegralFactory>& fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}

RelPotentialSOInt::RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>& aoint, const IntegralFactory* fact)
    : OneBodySOInt(aoint, fact) {
    natom_ = ob_->basis1()->molecule()->natom();
}
