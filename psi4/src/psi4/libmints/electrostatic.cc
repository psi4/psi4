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

#include "psi4/libmints/electrostatic.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

#include <libint2/engine.h>

using namespace psi;
using Zxyz_vector = std::vector<std::pair<double, std::array<double, 3>>>;

// Initialize potential_recur_ to +1 basis set angular momentum
ElectrostaticInt::ElectrostaticInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                   std::shared_ptr<BasisSet> bs2, int deriv)
    : PotentialInt(st, bs1, bs2, deriv) {
    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::nuclear, max_nprim, max_am, 0);

    if (deriv > 0) {
        throw PSIEXCEPTION("ElectrostaticInt: Derivatives are not supported");
    }
    set_chunks(1);

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

ElectrostaticInt::~ElectrostaticInt() {}

void ElectrostaticInt::compute(SharedMatrix& result, const Vector3& C) {
    engine0_->set_params(std::vector<std::pair<double, std::array<double,3>>>{{1.0, {C[0], C[1], C[2]}}});
    OneBodyAOInt::compute(result);
}

void ElectrostaticInt::set_origin(const Vector3& _origin) {
    origin_ = _origin;
    Zxyz_vector pcs;
    // l2 includes the electron charge, legacy psi4 code does not, so
    // we adopt this behavior here with -1.0
    pcs.push_back({-1.0, {origin_[0], origin_[1], origin_[2]}});
    engine0_->set_params(pcs);
}

SharedVector ElectrostaticInt::nuclear_contribution(std::shared_ptr<Molecule> mol) {
    auto sret = std::make_shared<Vector>(mol->natom());
    double* ret = sret->pointer();

    int natom = mol->natom();
    for (int k = 0; k < natom; k++) {
        Vector3 kgeom = mol->xyz(k);
        for (int i = 0; i < natom; i++) {
            if (i != k) {
                Vector3 igeom = mol->xyz(i);

                double x = kgeom[0] - igeom[0];
                double y = kgeom[1] - igeom[1];
                double z = kgeom[2] - igeom[2];
                double r2 = x * x + y * y + z * z;
                double r = std::sqrt(r2);
                ret[k] += mol->Z(i) / r;
            }
        }
    }

    return sret;
}
