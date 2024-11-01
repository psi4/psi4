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
    int current_ecp_number_;

    /// A list of all of the orbital basis sets, converted to LibECP format for the bra.
    std::vector<libecpint::GaussianShell> libecp_shells1_;
    /// A list of all of the orbital basis sets, converted to LibECP format for the ket.
    std::vector<libecpint::GaussianShell> libecp_shells2_;
    /// A list of the centers with ECPs, and the corresponding LibECP data structure.
    std::vector<std::pair<int,libecpint::ECP>> centers_and_libecp_ecps_;
    /// Tracks the iterations over ECP-bearing centers in Hessian integral calculations.
    int current_ecp_iterator_ = -1;
   public:
    ECPInt(std::vector<SphericalTransform> &, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv = 0);
    ~ECPInt() override;

    /// Initialize the iterator over ECPs for Hessian integral calculations (to save memory).  Call this first
    /// and then loop over ECP centers by calling `next_hessian()` in a while loop before calling the
    /// `compute_pair_deriv1(s1,s2)` function.
    void setup_hessian_iterations() { current_ecp_iterator_ = -1; }
    /// Advances to the next ECP center in Hessian integral calculations; this is designed to be used in a while
    /// loop and will return false when the iteration is finished.  Make sure that `setup_hessian_iterations()`
    /// is called before the while loop that calls this function.
    bool next_hessian_ecp() { return (++current_ecp_iterator_) != centers_and_libecp_ecps_.size(); }
    /// When looping over ECPs for Hessian integral calculations, call this function to find out which
    /// center hold the ECP corresponding to the current perturbation in the iterator.
    int current_ecp_center() const { return centers_and_libecp_ecps_[current_ecp_iterator_].first; }

    /// Overridden shell-pair integral calculation over all ECP centers
    void compute_shell(int s1, int s2) override;
    void compute_shell_deriv1(int s1, int s2) override;
    void compute_shell_deriv2(int s1, int s2) override;
};

class ECPSOInt : public OneBodySOInt {
    int natom_;

   public:
    ECPSOInt(const std::shared_ptr<OneBodyAOInt> &, const std::shared_ptr<IntegralFactory> &);
    ECPSOInt(const std::shared_ptr<OneBodyAOInt> &, const IntegralFactory *);
};

}  // namespace psi
#endif
