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

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/mcmurchiedavidson.h"

namespace psi {

class BasisSet;
class OneBodyAOInt;
class SphericalTransform;

/*! \ingroup MINTS
 *  \class MultipolePotentialInt
 *  \brief Computes multipole potential integrals, needed for EFP/PE calculations.
 *
 *  Use an IntegralFactory to create this object.
 *  The compute method takes a vector of SharedMatrix objects, which will be populated
 *  in CCA lexicographic (alphabetical) order of Cartesian components.
 */
class MultipolePotentialInt : public OneBodyAOInt, public mdintegrals::MDHelper {
    // maximum multipole potential order to compute (order of the 1/R derivative)
    int order_;

    //! CCA-ordered Cartesian components for the multipoles
    std::vector<std::vector<std::array<int, 3>>> comps_der_;

    //! Boys function evaluator from Libint2
    std::shared_ptr<const libint2::FmEval_Chebyshev7<double>> fm_eval_;

    //! R matrix (9.5.31)
    std::vector<double> R;

    //! Computes the multipole potential between two Gaussian shells.
    void compute_pair(const libint2::Shell&, const libint2::Shell&) override;

   public:
    //! Constructor. Do not call directly use an IntegralFactory.
    MultipolePotentialInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
                          int order, int deriv = 0);
    //! Virtual destructor
    ~MultipolePotentialInt() override;
};

}  // namespace psi
