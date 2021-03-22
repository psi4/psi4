/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2016-2017 Robert A. Shaw.
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

/* 	class ECPInt calculates one-body ECP integrals
        class AngularIntegral calculates and stores the angular integrals needed for the ECP integration
        class RadialIntegral abstracts the calculation of the radial integrals needed for the ECP integration,
        such that if a different approach was desired later, this could be done with minimal alterations to ECPInt.

        Robert A. Shaw 2016

        REFERENCES:
        (Flores06) R. Flores-Moreno et al., J. Comput. Chem. 27 (2006), 1009-1019
        (MM81) L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44 (1981), 289 - 301
*/

#ifndef LIBMINTS_ECPINT_H
#define LIBMINTS_ECPINT_H

#include <map>
#include <vector>

#include "psi4/pragma.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/sointegral_onebody.h"

#include "libecpint/ecpint.hpp"

namespace psi {

class BasisSet;
class GaussianShell;
class IntegralFactory;
class SphericalTransform;

/**
 * \ingroup MINTS
 * \class ECPIntegral
 * \brief Calculates ECP integrals.
 *
 * Given an ECP basis, and orbital bases, this will calculate the ECP integrals over all ECP centers.
 * The underlying engine is the libECPInt library (https://github.com/robashaw/libecpint).
 */
class ECPInt : public OneBodyAOInt {
   private:
    /// The LibECP instance that will do all the heavy lifting
    libecpint::ECPIntegral engine_;

    /// Integrals are requested using Psi4 GaussianShell objects, but we need the corresponding LibECP
    /// object to actually compute them.  This structure maps the two, via the Psi4 GaussianShell's
    /// first() function, which provides the index of the first basis function the shell involves.
    std::map<int, int> libecp_shell_lookup_;
    std::vector<libecpint::GaussianShell> libecp_shells_;
    std::vector<std::pair<int,libecpint::ECP>> centers_and_libecp_ecps_;

   public:
    ECPInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv = 0);
    ~ECPInt() override;
    void compute_pair(const libint2::Shell &shellA, const libint2::Shell &shellB) override;
    void compute_pair_deriv1(const libint2::Shell &shellA, const libint2::Shell &shellB) override;
};

class ECPSOInt : public OneBodySOInt {
    int natom_;

   public:
    ECPSOInt(const std::shared_ptr<OneBodyAOInt> &, const std::shared_ptr<IntegralFactory> &);
    ECPSOInt(const std::shared_ptr<OneBodyAOInt> &, const IntegralFactory *);
};

}  // namespace psi
#endif
