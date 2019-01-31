/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_electricfield_h_
#define _psi_src_lib_libmints_electricfield_h_

#include <vector>
#include "typedefs.h"
#include "psi4/libmints/osrecur.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/integral.h"
namespace psi {
class Molecule;

/*! \ingroup MINTS
 *  \class ElectricFieldInt
 *  \brief Computes electric field integrals.
 *
 *  Use an IntegralFactory to create this object.
 */
class ElectricFieldInt : public OneBodyAOInt {
    //! Obara and Saika recursion object to be used.
    ObaraSaikaTwoCenterElectricField efield_recur_;

    //! Number of atoms.
    int natom_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const GaussianShell &, const GaussianShell &) override;

    //! Computes the electric field gradient between two gaussian shells.
    void compute_pair_deriv1(const GaussianShell &, const GaussianShell &) override;

   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    ElectricFieldInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                     int deriv = 0);
    //! Virtual destructor
    ~ElectricFieldInt() override;

    //! Does the method provide first derivatives?
    bool has_deriv1() override { return false; }

    static Vector3 nuclear_contribution(const Vector3 &origin, std::shared_ptr<Molecule> mol);
    static SharedMatrix nuclear_contribution_to_gradient(const Vector3 &origin, std::shared_ptr<Molecule> mol);
};

}  // namespace psi

#endif
