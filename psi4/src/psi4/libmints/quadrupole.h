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

#ifndef _psi_src_lib_libmints_quadrupole_h_
#define _psi_src_lib_libmints_quadrupole_h_

#include <vector>
#include "typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/integral.h"
namespace psi {

class Molecule;

/*! \ingroup MINTS
 *  \class QuadrupoleInt
 *  \brief Computes quadrupole integrals. At last check this may not be working.
 *  Use an IntegralFactory to create this object.
 */
class QuadrupoleInt : public OneBodyAOInt {

   public:
    QuadrupoleInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>);
    ~QuadrupoleInt() override;

    void compute_pair(const libint2::Shell &, const libint2::Shell &) override;

    static SharedVector nuclear_contribution(std::shared_ptr<Molecule> mol, const Vector3 &origin);
};

}  // namespace psi

#endif
