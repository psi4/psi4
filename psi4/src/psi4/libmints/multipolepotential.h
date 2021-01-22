/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_multipolepotential_h_
#define _psi_src_lib_libmints_multipolepotential_h_

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/osrecur.h"

namespace psi {

class BasisSet;
class GaussianShell;
class OneBodyAOInt;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class MultipolePotentialInt
 *  \brief Computes multipole potential integrals, needed for EFP/PE calculations.
 *  Currently computes potential integrals up to order max_k,
 *  where the maximum is max_k = 3 (octupoles)
 *
 *  Use an IntegralFactory to create this object.
 *  The compute method takes a vector of SharedMatrix objects, which will be populated
 *  in alphabetical order of Cartesian components.
 */
class MultipolePotentialInt : public OneBodyAOInt {
    // OS Recursion for this type of potential integral
    ObaraSaikaTwoCenterMultipolePotentialRecursion mvi_recur_;

    // maximum multipole order to compute
    int max_k_;

    //! Computes the electric field between two gaussian shells.
    void compute_pair(const GaussianShell&, const GaussianShell&) override;

   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    MultipolePotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                          int max_k = 0, int deriv = 0);
    //! Virtual destructor
    ~MultipolePotentialInt() override;
};

}  // namespace psi

#endif
