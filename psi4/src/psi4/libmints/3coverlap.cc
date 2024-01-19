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
#include "psi4/libmints/3coverlap.h"
#include "psi4/libmints/basisset.h"

#include <memory>
#include <stdexcept>

#include <libint2/engine.h>

using namespace psi;

ThreeCenterOverlapInt::ThreeCenterOverlapInt(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2,
                                             std::shared_ptr<BasisSet> bs3)
    : bs1_(bs1), bs2_(bs2), bs3_(bs3) {
    int maxam1 = bs1_->max_am();
    int maxam2 = bs2_->max_am();
    int maxam3 = bs3_->max_am();
    int max_am = std::max({maxam1, maxam2, maxam3});
    int max_nprim = std::max({basis1()->max_nprimitive(), basis2()->max_nprimitive(), basis3()->max_nprimitive()});

    int max_nao = INT_NCART(maxam1) * INT_NCART(maxam2) * INT_NCART(maxam3);
    zero_vec_ = std::vector<double>(max_nao, 0.0);

    // set engine precision to 0.0 to disable primitive screening
    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::delta, max_nprim, max_am, 0, 0.0);
    buffers_.resize(1);
}

ThreeCenterOverlapInt::~ThreeCenterOverlapInt() {}

std::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis1() { return bs1_; }

std::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis2() { return bs2_; }

std::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis3() { return bs3_; }

void ThreeCenterOverlapInt::compute_shell(int sh1, int sh2, int sh3) {
    compute_pair(bs1_->l2_shell(sh1), bs2_->l2_shell(sh2), bs3_->l2_shell(sh3));
}

void ThreeCenterOverlapInt::compute_pair(const libint2::Shell& s1, const libint2::Shell& s2, const libint2::Shell& s3) {
    engine0_->compute(s1, s2, s3, libint2::Shell::unit());
    buffers_[0] = engine0_->results()[0];
    // in case l2 gives us a nullptr, point to a zero vector
    // to avoid undefined behavior in the calling code
    if (buffers_[0] == nullptr) {
        buffers_[0] = zero_vec_.data();
    }
}
