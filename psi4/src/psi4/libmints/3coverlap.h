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
#pragma once

#include "psi4/libmints/integral.h"

namespace libint2 {
class Engine;
class Shell;
}  // namespace libint2

namespace psi {

class BasisSet;

/** \ingroup MINTS
    \class ThreeCenterOverlapInt
    \brief Three center overlap integral.
 */
class ThreeCenterOverlapInt {
   protected:
    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::shared_ptr<BasisSet> bs3_;

    std::unique_ptr<libint2::Engine> engine0_;

    /// Buffer to hold the source integrals.
    std::vector<const double*> buffers_;

    /// A vector of zeros that we can point to if libint2 gives us back a nullptr
    std::vector<double> zero_vec_;

    void compute_pair(const libint2::Shell& s1, const libint2::Shell& s2, const libint2::Shell& s3);

   public:
    ThreeCenterOverlapInt(std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, std::shared_ptr<BasisSet> bs3);
    ~ThreeCenterOverlapInt();

    /// Basis set on center one.
    std::shared_ptr<BasisSet> basis1();
    /// Basis set on center two.
    std::shared_ptr<BasisSet> basis2();
    /// Basis set on center three.
    std::shared_ptr<BasisSet> basis3();

    /// Compute the integrals of the form (a|c|b).
    virtual void compute_shell(int, int, int);

    /// Buffer where each chunk of integrals is placed
    const std::vector<const double*>& buffers() const { return buffers_; }
};

}  // namespace psi