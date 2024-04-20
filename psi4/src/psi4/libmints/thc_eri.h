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

#ifndef PSI4_SRC_LIBMINTS_LS_THC_H_
#define PSI4_SRC_LIBMINTS_LS_THC_H_

#include "basisset.h"
#include "matrix.h"

#include <vector>

namespace psi {

class THC_Computer {
   protected:
    // The molecule to evaluate the THC factors
    std::shared_ptr<Molecule> molecule_;
    // The primary basis set
    std::shared_ptr<BasisSet> primary_;

    // THC Factor for first index
    SharedMatrix x1_;
    // THC Factor for second index
    SharedMatrix x2_;
    // THC Factor for third index
    SharedMatrix x3_;
    // THC Factor for fourth index
    SharedMatrix x4_;
    // THC connecting factor
    SharedMatrix Z_IJ_;

   public:
    THC_Computer(std::shared_ptr<Molecule> molecule);
    virtual ~THC_Computer();

    /// Compute THC Factors
    void compute_thc_factorization();

    /// Returns the molecule the THC factors are evaluated on
    std::shared_ptr<Molecule> molecule() const { return molecule_; }
    /// Returns the THC factor for the first index
    SharedMatrix get_x1() const { return x1_; }
    /// Returns the THC factor for the second index
    SharedMatrix get_x2() const { return x2_; }
    /// Returns the THC factor for the third index
    SharedMatrix get_x3() const { return x3_; }
    /// Returns the THC factor for the fourth index
    SharedMatrix get_x4() const { return x4_; }

};

// Least Squares Tensor Hypercontraction
class LS_THC_Computer : public THC_Computer {
   protected:
    /// Parrish LS-THC Algorithm X
    void build_E_df();
    /// Parrish LS-THC Algorithm Y
    void build_E_exact();

   public:
    LS_THC_Computer(std::shared_ptr<Molecule> molecule);
    ~LS_THC_Computer() override;
};

}