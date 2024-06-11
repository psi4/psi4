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

#include <vector>
#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/libmints/onebody.h"

namespace psi {

class BasisSet;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class PotentialErfInt
 *  \brief Computes the erf-attenuated Coulomb integrals on a given origin.
 * Use an IntegralFactory to create this object.
 */
class PotentialErfInt : public OneBodyAOInt {
   protected:
    /// The range-separation parameter. Defaults to 1.0
    double omega_;

   public:
    /// Constructor
    PotentialErfInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                    double omega = 1.0, int deriv = 0);
    ~PotentialErfInt() override{};

    void set_origin(const Vector3& _origin) override;
};

/*! \ingroup MINTS
 *  \class PotentialErfComplementInt
 *  \brief Computes the complementary erf-attenuated Coulomb integrals on a given origin.
 * Use an IntegralFactory to create this object.
 */

class PotentialErfComplementInt : public OneBodyAOInt {
   protected:
    /// The range-separation parameter. Defaults to 1.0
    double omega_;

   public:
    /// Constructor
    PotentialErfComplementInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                              double omega = 1.0, int deriv = 0);
    ~PotentialErfComplementInt() override{};

    void set_origin(const Vector3& _origin) override;
};

}  // namespace psi
