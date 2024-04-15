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

#ifndef _psi_src_lib_libmints_rel_potential_h_
#define _psi_src_lib_libmints_rel_potential_h_

#include <vector>
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/sointegral_onebody.h"

namespace psi {
// TODO:  This is all in typedefs.h ....
class BasisSet;
class OneBodyAOInt;
class IntegralFactory;
class SphericalTransform;
class OneBodySOInt;
class CdSalcList;

/*! \ingroup MINTS
 *  \class RelPotentialInt
 *  \brief Computes relativistic potential integrals.
 * Use an IntegralFactory to create this object.
 */
class RelPotentialInt : public OneBodyAOInt {

   protected:
    /// The charges and locations that define the external potential
    std::vector<std::pair<double, std::array<double, 3>>> Zxyz_;

   public:
    /// Constructor. Assumes nuclear centers/charges as the potential
    RelPotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                    int deriv = 0);

    /// Set the field of charges
    void set_charge_field(const std::vector<std::pair<double, std::array<double, 3>>>& Zxyz);
};

class RelPotentialSOInt : public OneBodySOInt {
    int natom_;

   public:
    RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>&, const std::shared_ptr<IntegralFactory>&);
    RelPotentialSOInt(const std::shared_ptr<OneBodyAOInt>&, const IntegralFactory*);
};

}  // namespace psi

#endif
