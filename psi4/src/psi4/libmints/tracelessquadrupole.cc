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
#include "psi4/libmints/tracelessquadrupole.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libciomr/libciomr.h"

#include <libint2/engine.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

using namespace psi;

// Initialize overlap_recur_ to +2 basis set angular momentum
TracelessQuadrupoleInt::TracelessQuadrupoleInt(std::vector<SphericalTransform>& st, std::shared_ptr<BasisSet> bs1,
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

TracelessQuadrupoleInt::~TracelessQuadrupoleInt() { }

void TracelessQuadrupoleInt::compute_pair(const libint2::Shell &s1, const libint2::Shell &s2) {
    engine0_->compute(s1, s2);

    size_t nints = s1.size() * s2.size();

    // Computed ints are now
    //
    // 0  1  2  3   4   5   6   7   8   9
    // s µx µy µz  Qxx Qxy Qxz Qyy Qyz Qzz
    //
    for (int i = 0; i < nints; ++i) {
        double Qxx = -engine0_->results()[4][i];
        double Qxy = -engine0_->results()[5][i];
        double Qxz = -engine0_->results()[6][i];
        double Qyy = -engine0_->results()[7][i];
        double Qyz = -engine0_->results()[8][i];
        double Qzz = -engine0_->results()[9][i];
        double R2 = (Qxx + Qyy + Qzz) / 3;
        const_cast<double*>(engine0_->results()[0])[i] = 1.5 * (Qxx - R2);
        const_cast<double*>(engine0_->results()[1])[i] = 1.5 * Qxy;
        const_cast<double*>(engine0_->results()[2])[i] = 1.5 * Qxz;
        const_cast<double*>(engine0_->results()[3])[i] = 1.5 * (Qyy - R2);
        const_cast<double*>(engine0_->results()[4])[i] = 1.5 * Qyz;
        const_cast<double*>(engine0_->results()[5])[i] = 1.5 * (Qzz - R2);
    }

    for (int chunk = 0; chunk < 6; ++chunk) {
        buffers_[chunk] = engine0_->results()[chunk];
    }
}
