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

#ifndef _psi_src_lib_libmints_overlap_h_
#define _psi_src_lib_libmints_overlap_h_

#include <vector>
#include "psi4/libmints/onebody.h" // needed because we derive from OneBodyAOInt
#include "psi4/libmints/osrecur.h"

namespace psi {

    class BasisSet;
    class GaussianShell;
    class IntegralFactory;
    class SphericalTransform;
    class Matrix;

/*! \ingroup MINTS
 *  \class OverlapInt
 *  \brief This class computes overlap integrals and soon overlap integral derivatives.
 *  Use an IntegralFactory to create this object.
 */
class OverlapInt : public OneBodyAOInt
{
    /// Generic Obara Saika recursion object.
    ObaraSaikaTwoCenterRecursion overlap_recur_;

    /// Computes the overlap between a given shell pair.
    void compute_pair(const GaussianShell& , const GaussianShell&);
    void compute_pair_deriv1(const GaussianShell& s1, const GaussianShell& s2);
    void compute_pair_deriv2(const GaussianShell&, const GaussianShell&);

public:
    /// Constructor, it assumes you are not computing derivatives by default
    OverlapInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    virtual ~OverlapInt();

    /// Does the method provide first derivatives?
    bool has_deriv1() { return true; }
    /// Does the method provide second derivatives?
    bool has_deriv2() { return true; }
};

}

#endif
