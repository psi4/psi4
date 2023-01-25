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
#include <algorithm>
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

#include <libint2/engine.h>

using namespace psi;

DipoleInt::DipoleInt(std::vector<SphericalTransform> &spherical_transforms, std::shared_ptr<BasisSet> bs1,
                     std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv) {
    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    if (nderiv == 0) {
        set_chunks(3);
        engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::emultipole1, max_nprim, max_am, 0);
    } else {
        throw PSIEXCEPTION("No derivatives available from DipoleInt. Use MultipoleInt instead.");
    }
    buffers_.resize(nchunk_);
}

DipoleInt::~DipoleInt() {}

SharedVector DipoleInt::nuclear_contribution(std::shared_ptr<Molecule> mol, const Vector3 &origin) {
    auto sret = std::make_shared<Vector>(3);
    double *ret = sret->pointer();

    for (int i = 0; i < mol->natom(); ++i) {
        Vector3 geom = mol->xyz(i) - origin;
        ret[0] += mol->Z(i) * geom[0];
        ret[1] += mol->Z(i) * geom[1];
        ret[2] += mol->Z(i) * geom[2];
    }

    return sret;
}

SharedMatrix DipoleInt::nuclear_gradient_contribution(std::shared_ptr<Molecule> mol) {
    auto sret = std::make_shared<Matrix>("Nuclear dipole derivative (3Nx3)", 3 * mol->natom(), 3);
    double **ret = sret->pointer();

    for (int i = 0; i < mol->natom(); ++i) {
        ret[3 * i + 0][0] = mol->Z(i);
        ret[3 * i + 1][1] = mol->Z(i);
        ret[3 * i + 2][2] = mol->Z(i);
    }

    return sret;
}

void DipoleInt::compute_pair(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine0_->compute(s1, s2);

    size_t nints = s1.size() * s2.size();
    // Libint gives us the overlap, mu_x, mu_y, mu_z in the buffers.
    // We don't care about the overlap here so we just skip over it.
    for (int chunk = 1; chunk < 4; chunk++) {
        double *ptr = const_cast<double *>(engine0_->results()[chunk]);
        std::transform(ptr, ptr + nints, ptr, [](double val) -> double { return -val; });
        buffers_[chunk - 1] = engine0_->results()[chunk];
    }
}