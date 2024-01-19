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

#ifndef _psi_src_lib_libmints_kinetic_h_
#define _psi_src_lib_libmints_kinetic_h_

#include "psi4/pragma.h"
#include <memory>
#include <vector>
#include "psi4/libmints/onebody.h"

namespace psi {

class BasisSet;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class KineticInt
 *  \brief Computes kinetic integrals.
 *
 * Use an IntegralFactory to create this object.
 */
class KineticInt : public OneBodyAOInt {

   public:
    //! Constructor. Do not call directly, use an IntegralFactory.
    KineticInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv = 0);
    //! Virtual destructor.
    ~KineticInt() override;

    /// Does the method provide first derivatives?
    bool has_deriv1() override { return true; }

    /// Does the method provide first derivatives?
    bool has_deriv2() override { return true; }

};

}  // namespace psi

#endif
