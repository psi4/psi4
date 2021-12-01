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
#include "psi4/libmints/dipole.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libciomr/libciomr.h"

#include <libint2/engine.h>

using namespace psi;

// Initialize overlap_recur_ to +1 basis set angular momentum, +1 on each center is sufficient
// to compute the dipole derivatives
DipoleInt::DipoleInt(std::vector<SphericalTransform> &spherical_transforms, std::shared_ptr<BasisSet> bs1,
                     std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv) {
    // int maxam1 = bs1_->max_am();
    // int maxam2 = bs2_->max_am();

    // int maxnao1 = (maxam1 + 1) * (maxam1 + 2) / 2;
    // int maxnao2 = (maxam2 + 1) * (maxam2 + 2) / 2;

    // // Increase buffer size to handle x, y, and z components
    // if (deriv_ == 0) {
    //     buffer_ = new double[3 * maxnao1 * maxnao2];
    //     set_chunks(3);
    // } else if (deriv_ == 1) {
    //     natom_ = bs1_->molecule()->natom();
    //     buffer_ = new double[6 * 3 * maxnao1 * maxnao2];
    //     set_chunks(6 * 3);
    // }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    if (nderiv == 0) {
        set_chunks(3);

        engine0_ =
            std::unique_ptr<libint2::Engine>(new libint2::Engine(libint2::Operator::emultipole1, max_nprim, max_am, 0));
    } else if (nderiv == 1) {
        // We set chunk count for normalize_am and pure_transform
        set_chunks(18);

        engine0_ =
            std::unique_ptr<libint2::Engine>(new libint2::Engine(libint2::Operator::emultipole1, max_nprim, max_am, 0));
        engine1_ =
            std::unique_ptr<libint2::Engine>(new libint2::Engine(libint2::Operator::emultipole1, max_nprim, max_am, 1));
    }

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

DipoleInt::~DipoleInt() { delete[] buffer_; }

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

    for (int chunk = 1; chunk < 4; chunk++) {
        buffers_[chunk - 1] = engine0_->results()[chunk];
    }
}

void DipoleInt::compute_pair_deriv1(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine1_->compute(s1, s2);

    for (int i = 0; i < 6; i) {
        buffers_[3 * i + 0] = engine1_->results()[4 * i + 1];
        buffers_[3 * i + 1] = engine1_->results()[4 * i + 2];
        buffers_[3 * i + 2] = engine1_->results()[4 * i + 3];
    }
}
