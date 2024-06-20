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
#include "molecule.h"
#include "matrix.h"
#include "psi4/liboptions/liboptions.h"

#include <vector>

namespace psi {

class THC_Computer {
   protected:
    // The molecule to evaluate the THC factors
    std::shared_ptr<Molecule> molecule_;
    // The primary basis set
    std::shared_ptr<BasisSet> primary_;
    // Options object
    Options& options_;

    // THC Factor for first index
    SharedMatrix x1_;
    // THC Factor for second index
    SharedMatrix x2_;
    // THC Factor for third index
    SharedMatrix x3_;
    // THC Factor for fourth index
    SharedMatrix x4_;
    // THC connecting factor
    SharedMatrix Z_PQ_;

   public:
    THC_Computer(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, Options& options);
    virtual ~THC_Computer();

    /// Compute THC Factors
    virtual void compute_thc_factorization() = 0;

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
    /// Returns the THC connecting factor
    SharedMatrix get_Z() const { return Z_PQ_; }

};

// Least Squares Tensor Hypercontraction
// Derived from work of Parrish et al. 2012 (doi: 10.1063/1.4768233)
class LS_THC_Computer : public THC_Computer {
   protected:
    /// Use DF integrals to perform LS-THC ?
    bool use_df_;
    
    /// Auxiliary basis set (null if not using DF approximation)
    std::shared_ptr<BasisSet> auxiliary_;

    /// Print options and other info for LS-THC decomposition
    void print_header();

    /// Parrish LS-THC Procedure 2
    SharedMatrix build_E_exact();
    /// Parrish LS-THC Procedure 3
    SharedMatrix build_E_df();
    /// Matthews 2020 SI Page 5 (returns S matrix using rank-reduced grid)
    /// doi: 10.1021/acs.jctc.9b01205
    SharedMatrix prune_grid();

   public:
    LS_THC_Computer(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, Options& options);
    LS_THC_Computer(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    ~LS_THC_Computer() override;

    /// Compute THC Factors using LS-THC factorization
    void compute_thc_factorization() override;
};

}

#endif