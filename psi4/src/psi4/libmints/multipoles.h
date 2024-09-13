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
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/mcmurchiedavidson.h"

namespace psi {
class Molecule;

/*! \ingroup MINTS
 *  \class MultipoleInt
 *  \brief Computes arbitrary-order multipole integrals.
 *
 * Use an IntegralFactory to create this object. */
class MultipoleInt : public OneBodyAOInt, public mdintegrals::MDHelper {
    //! The order of multipole moment to compute
    int order_;

    //! Multipole intermediates
    //! M matrix (9.5.31)
    std::vector<double> Mx;
    std::vector<double> My;
    std::vector<double> Mz;
    //! S matrix (9.5.29)
    std::vector<double> Sx;
    std::vector<double> Sy;
    std::vector<double> Sz;

    //! CCA-ordered Cartesian components for the multipoles
    std::vector<std::vector<std::array<int, 3>>> comps_mul_;

   public:
    //! Constructor. Do not call directly. Use an IntegralFactory.
    MultipoleInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int order,
                 int deriv = 0);
    //! Virtual destructor
    ~MultipoleInt() override;

    //! Computes the multipole integrals between two gaussian shells.
    void compute_pair(const libint2::Shell &, const libint2::Shell &) override;

    //! Computes the first derivative of the multipole integrals between two gaussian shells.
    void compute_pair_deriv1(const libint2::Shell &, const libint2::Shell &) override;

    //! Does the method provide first derivatives?
    bool has_deriv1() override { return true; }

    /// Returns the nuclear contribution to the multipole moments, with angular momentum up to order
    static SharedVector nuclear_contribution(std::shared_ptr<Molecule> mol, int order, const Vector3 &origin);
};

}  // namespace psi