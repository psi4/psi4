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

#ifndef _psi_src_lib_libmints_pseudospectral_h_
#define _psi_src_lib_libmints_pseudospectral_h_

#include <vector>
 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/osrecur.h"

namespace psi {

    class BasisSet;
    class GaussianShell;
    class IntegralFactory;
    class SphericalTransform;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes pseudospectral integrals.
 * Use an IntegralFactory to create this object.
 */
class PseudospectralInt : public OneBodyAOInt
{
    /// Computes integrals between two shell objects.
    void compute_pair(const GaussianShell&, const GaussianShell&);
    /// Computes integrals between two shell objects.
    void compute_pair_deriv1(const GaussianShell&, const GaussianShell&);

protected:

    /// Use range-separation or not? Defaults to false. If so, produce <m|erf(\omega r) / r|n> integrals
    bool use_omega_;

    /// The range-separation parameter. Defaults to 0.0
    double omega_;

    /// The integration point
    double C_[3];
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIRecursion potential_recur_;
    /// Recursion object that does the heavy lifting.
    ObaraSaikaTwoCenterVIDerivRecursion potential_deriv_recur_;

public:
    /// Constructor
    PseudospectralInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    ~PseudospectralInt();

    /// Computes integrals between two shells.
    void compute_shell_deriv1(int, int);

    /// Set integration point
    void set_point(double x, double y, double z) { C_[0] = x; C_[1] = y; C_[2] = z; }

    /// Set omega value, turns use_omega_ to true
    void set_omega(double omega) { use_omega_ = (omega != 0.0); omega_ = omega; }

    /// Set the value of the use_omega_ flag
    void use_omega(bool yes) { use_omega_ = yes; }

    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }
};

}

#endif
