/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#ifndef SPLITJK_H
#define SPLITJK_H

#include <vector>

#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
#include <unordered_map>
PRAGMA_WARNING_POP
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/dimension.h"
//#include "jk.h"

namespace psi {
class MinimalInterface;
class BasisSet;
class Matrix;
class ERISieve;
class TwoBodyAOInt;
class Options;
class PSIO;
class DFHelper;
class DFTGrid;

namespace pk {
class PKManager;
}

/**
 * Class SplitJK 
 *
 * The base class defining the backend to CompositeJK. 
 * Each derived class for SplitJK contains one of the
 * subalgorithms that can be used within the SplitJK
 * framework. 
 * 
 * Current algorithms in place:
 * J: Direct DF-J
 * K: COSX, LinK
 *
 */
class PSI_API SplitJK : public JK {
   protected:

    /// The number of threads to be used for integral computation
    int nthreads_;
    /// Options object
    Options& options_;

    /// SplitJK algorithm info
    std::string algo_; 

    // Perform Density matrix-based integral screening?
    bool density_screening_;

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    SplitJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    /// Destructor
    ~SplitJK() override;

    /// Build either the coulomb (J) matrix or the exchange (K) matrix
    /// using a given algorithm 
    void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp) = 0;

    // => Knobs <= //
    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() = 0; 
 
    /**
    * print name of method
    */ 
    std::string name() = 0; 
}

// ==> Start SplitJK Coulomb (J) Algorithms here <==

/// Build the coulomb (J) matrix using Direct DF-J
/// Reference is https://doi.org/10.1039/B204199P
class PSI_API DirectDFJ : public SplitJK {
    // => Density Fitting Stuff <= //

    /// Auxiliary basis set
    std::shared_ptr<BasisSet> auxiliary_;
    /// Coulomb Metric
    SharedMatrix J_metric_;
    /// per-thread TwoBodyAOInt object (for computing three/four-center ERIs)
    std::unordered_map<std::string, std::vector<std::shared_ptr<TwoBodyAOInt>>> eri_computers_;

    /// Common initialization
    void common_init();

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    DirectDFJ(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    /// Destructor
    ~DirectDFJ() override;

   void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp) override;

    // => Knobs <= //
    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() const override;

    size_t num_computed_shells() {
        outfile->Printf("WARNING: JK::num_computed_shells() was called, but benchmarking is disabled for the chosen JK algorithm.");
        outfile->Printf(" Returning 0 as computed shells count.\n");

        return 0;
    }

    /**
    * print name of method
    */ 
    std::string name() override { return "DirectDFJ"; }
};

// ==> Start SplitJK Exchange (K) Algorithms here <==

/**
 * @author Andy Jiang, Georgia Tech, December 2021
 * 
 * @brief constructs the K matrix using the LinK algorithm, described in [Ochsenfeld:1998:1663]_
 * doi: 10.1063/1.476741
 * 
 * @param ints A list of TwoBodyAOInt objects (one per thread) to optimize parallel efficiency
 * @param D The list of AO density matrices to contract to form J and K (1 for RHF, 2 for UHF/ROHF)
 * @param K The list of AO K matrices to build (Same size as D)
 * 
 */
#if 0
class PSI_API LinK : public SplitJK {
    // => LinK variables <= //

    // Density-based ERI Screening tolerance to use in the LinK algorithm
    double linK_ints_cutoff_;

    /// Common initialization
    void common_init();

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    LinK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    /// Destructor
    ~LinK() override;

    /// Build the exchange (K) matrix using LinK 
    void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp) override;

    // => Knobs <= //
    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() const override;
};

class PSI_API COSK : public SplitJK {
    // => Semi-Numerical Stuff <= //

    /// Small DFTGrid for initial SCF iterations
    std::shared_ptr<DFTGrid> grid_init_;
    /// Large DFTGrid for the final SCF iteration
    std::shared_ptr<DFTGrid> grid_final_;
    /// Overlap fitting metric for grid_initial_
    SharedMatrix Q_init_;
    /// Overlap fitting metric for grid_final_
    SharedMatrix Q_final_;
 
    /// Common initialization
    void common_init();

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    COSK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options);
    /// Destructor
    ~COSK() override;

    /// Build the exchange (K) matrix using COSX
    void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp) override;

    // => Knobs <= //
    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() const override;
};
#endif
}

#endif

