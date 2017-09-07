/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_electrostatic_h_
#define _psi_src_lib_libmints_electrostatic_h_

#include "psi4/libmints/potential.h"
#include "psi4/pragma.h"

namespace psi {

class BasisSet;
class Molecule;
class GaussianShell;
class SphericalTransform;
class Vector3;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class ElectrostaticInt : public PotentialInt
{
    void compute_pair(const GaussianShell&, const GaussianShell&)
    {}

public:
    /// Constructor
    ElectrostaticInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    ~ElectrostaticInt();

    // Intel C++ 12 thinks we're trying to overload the "void compute_shell(int, int)" and warns us about it.
    // The following line is to shut it up.
    #pragma warning disable 1125
    /// Computes integrals between two shells.
    void compute_shell(int, int, const Vector3&);
    /// Computes integrals between two shells.
    void compute_pair(const GaussianShell&, const GaussianShell&, const Vector3&);

    PRAGMA_WARNING_PUSH
    PRAGMA_WARNING_IGNORE_OVERLOADED_VIRTUAL
    /// Computes integrals and stores in result.
    void compute(SharedMatrix& result, const Vector3&);
    PRAGMA_WARNING_POP

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }

    static SharedVector nuclear_contribution(std::shared_ptr<Molecule> mol);
};

}

#endif
