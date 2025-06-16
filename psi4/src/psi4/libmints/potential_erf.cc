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
#include "psi4/libmints/potential_erf.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include <libint2/engine.h>

using namespace psi;
using Zxyz_vector = std::vector<std::pair<double, std::array<double, 3>>>;

PotentialErfInt::PotentialErfInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                 std::shared_ptr<BasisSet> bs2, double omega, int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv), omega_(omega) {
    if (deriv > 0) {
        throw PSIEXCEPTION("PotentialErfInt does not support derivatives.");
    }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::erf_nuclear, max_nprim, max_am, 0);

    // set engine parameters
    set_origin(Vector3(0, 0, 0));

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

void PotentialErfInt::set_origin(const Vector3& _origin) {
    origin_ = _origin;
    std::cout << "(set origin) erf using omega = " << omega_ << std::endl;
    Zxyz_vector pcs;
    // l2 includes the electron charge, legacy psi4 code does not, so
    // we adopt this behavior here with -1.0
    pcs.push_back({-1.0, {origin_[0], origin_[1], origin_[2]}});
    std::tuple<double, Zxyz_vector> params{omega_, pcs};
    engine0_->set_params(params);
}

PotentialErfComplementInt::PotentialErfComplementInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
                                                     std::shared_ptr<BasisSet> bs2, double omega, int deriv)
    : OneBodyAOInt(st, bs1, bs2, deriv), omega_(omega) {
    if (deriv > 0) {
        throw PSIEXCEPTION("PotentialErfComplementInt does not support derivatives.");
    }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::erfc_nuclear, max_nprim, max_am, 0);

    // set engine parameters
    set_origin(Vector3(0, 0, 0));

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

void PotentialErfComplementInt::set_origin(const Vector3& _origin) {
    origin_ = _origin;
    std::cout << "(set origin) erfc using omega = " << omega_ << std::endl;
    Zxyz_vector pcs;
    // l2 includes the electron charge, legacy psi4 code does not, so
    // we adopt this behavior here with -1.0
    pcs.push_back({-1.0, {origin_[0], origin_[1], origin_[2]}});
    std::tuple<double, Zxyz_vector> params{omega_, pcs};
    engine0_->set_params(params);
}
