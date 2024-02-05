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
#include "psi4/libmints/quadrupole.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libciomr/libciomr.h"

#include "libint2/engine.h"
;
using namespace psi;

// Initialize overlap_recur_ to +2 basis set angular momentum
QuadrupoleInt::QuadrupoleInt(std::vector<SphericalTransform> &st, std::shared_ptr<BasisSet> bs1,
                             std::shared_ptr<BasisSet> bs2)
    : OneBodyAOInt(st, bs1, bs2) {

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    if (deriv_ == 0) {
        set_chunks(6);

        engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::emultipole2, max_nprim, max_am, 0);
    } else {
        throw PSIEXCEPTION("Derivatives for quadrupole integrals are not implemented.");
    }

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

QuadrupoleInt::~QuadrupoleInt() { }

SharedVector QuadrupoleInt::nuclear_contribution(std::shared_ptr<Molecule> mol, const Vector3 &origin) {
    auto sret = std::make_shared<Vector>(6);
    double *ret = sret->pointer();

    for (int i = 0; i < mol->natom(); ++i) {
        Vector3 geom = mol->xyz(i) - origin;
        ret[0] += mol->Z(i) * geom[0] * geom[0];  // xx
        ret[1] += mol->Z(i) * geom[0] * geom[1];  // xy
        ret[2] += mol->Z(i) * geom[0] * geom[2];  // xz
        ret[3] += mol->Z(i) * geom[1] * geom[1];  // yy
        ret[4] += mol->Z(i) * geom[1] * geom[2];  // yz
        ret[5] += mol->Z(i) * geom[2] * geom[2];  // zz
    }

    return sret;
}

void QuadrupoleInt::compute_pair(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine0_->compute(s1, s2);

    size_t nints = s1.size() * s2.size();
    for (int chunk = 4; chunk < 10; chunk++) {
        double * ptr = const_cast<double*>(engine0_->results()[chunk]);
        std::transform(ptr, ptr + nints, ptr, [](double val) -> double { return -val; });
        buffers_[chunk - 4] = engine0_->results()[chunk];
    }
}
