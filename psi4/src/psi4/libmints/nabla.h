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

#ifndef _psi_src_lib_libmints_nabla_h_
#define _psi_src_lib_libmints_nabla_h_

#include "psi4/pragma.h"
#include <memory>
#include "psi4/libmints/onebody.h"

namespace psi {

class BasisSet;
class GaussianShell;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class NablaInt
 *  \brief Computes nabla integrals.
 *
 * Use an IntegralFactory to create this object. */
class NablaInt : public OneBodyAOInt {

   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    NablaInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv = 0);
    //! Virtual destructor
    ~NablaInt() override;

    bool is_antisymmetric() const override { return true; }
};

}  // namespace psi

#endif
