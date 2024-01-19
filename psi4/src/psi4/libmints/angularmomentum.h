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

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/mcmurchiedavidson.h"

namespace psi {

class BasisSet;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class AngularMomentumInt
 *  \brief Computes angular momentum integrals.
 *
 * Use an IntegralFactory to create this object. */
class AngularMomentumInt : public OneBodyAOInt, public mdintegrals::MDHelper {
   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    AngularMomentumInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                       int deriv = 0);
    //! Virtual destructor
    ~AngularMomentumInt() override;

    void compute_pair(const libint2::Shell &, const libint2::Shell &) override;

    bool is_antisymmetric() const override { return true; }
};

}  // namespace psi