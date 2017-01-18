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
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _psi_src_lib_libmints_efpmultipolepotential_h_
#define _psi_src_lib_libmints_efpmultipolepotential_h_

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/osrecur.h"

namespace psi {

class BasisSet;
class GaussianShell;
class ObaraSaikaTwoCenterRecursion;
class OneBodyAOInt;
class PotentialInt;
class IntegralFactory;
class SphericalTransform;
class Vector3;


/*! \ingroup MINTS
 *  \class MultipolePotentialInt
 *  \brief Computes multipole potential integrals, needed for EFP calculations.
 *  Currently computes potential integrals through octopoles.
 *
 *  Use an IntegralFactory to create this object.
 *  The compute method takes a vector of SharedMatrix objects, which will be populated
 *  in the following order, which matches the order expected by libEFP. The scale factors
 *  are the prefactors that should be included in any contractions involving these integrals
 *  where the numerator comes from permutational symmetry, and the denominator arises from
 *  the coefficients of the power series expansion of the Coulomb operator.

 *  Matrix # | Multipole Type | Scale Factor
 *  -----------------------------------------
 *           | // Charge      |
 *      0    |      0         |    1
 *           | // Dipole      |
 *      1    |      X         |    1
 *      2    |      Y         |    1
 *      3    |      Z         |    1
 *           | // Quadrupole  |
 *      4    |      XX        |    1/3
 *      5    |      YY        |    1/3
 *      6    |      ZZ        |    1/3
 *      7    |      XY        |    2/3
 *      8    |      XZ        |    2/3
 *      9    |      YZ        |    2/3
 *           | // Octupole    |
 *     10    |      XXX       |    1/15
 *     11    |      YYY       |    1/15
 *     12    |      ZZZ       |    1/15
 *     13    |      XXY       |    3/15
 *     14    |      XXZ       |    3/15
 *     15    |      XYY       |    3/15
 *     16    |      YYZ       |    3/15
 *     17    |      XZZ       |    3/15
 *     18    |      YZZ       |    3/15
 *     19    |      XYZ       |    6/15
 *
 */
class EFPMultipolePotentialInt : public OneBodyAOInt
{
    // OS Recursion for this type of potential integral
    ObaraSaikaTwoCenterEFPRecursion mvi_recur_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const GaussianShell&, const GaussianShell&);

public:
    //! Constructor. Do not call directly use an IntegralFactory.
    EFPMultipolePotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    //! Virtual destructor
    virtual ~EFPMultipolePotentialInt();

};

}

#endif
