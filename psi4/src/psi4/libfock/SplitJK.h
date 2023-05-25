/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

namespace psi {
class MinimalInterface;
class BasisSet;
class Matrix;
class TwoBodyAOInt;
class Options;
class PsiOutStream;
class DFTGrid;

/**
 * Class SplitJK
 *
 * The base class defining the backend to CompositeJK.
 * Each derived class for SplitJK contains one of the
 * subalgorithms that can be used within the CompositeJK
 * framework.
 *
 * Current algorithms in place:
 * J: DF-DirJ
 * K: COSX, LinK
 *
 */
class PSI_API SplitJK {
   protected:

    /// The number of threads to be used for integral computation
    int nthreads_;
    /// Options object
    Options& options_;

    /// Primary basis set
    std::shared_ptr<BasisSet> primary_;

    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    bool bench_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;

    // Perform Density matrix-based integral screening?
    bool density_screening_;
    /// Use severe screening techniques? Useful in early SCF iterations
    bool early_screening_;
    /// Left-right symmetric?
    bool lr_symmetric_;

    /// Number of ERI shell quartets computed, i.e., not screened out
    size_t num_computed_shells_;

   public:
    // => Constructors < = //
    SplitJK(std::shared_ptr<BasisSet> primary, Options& options);

    /// Destructor
    virtual ~SplitJK();

    /// Build either the coulomb (J) matrix or the exchange (K) matrix
    /// using a given algorithm
    virtual void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp,
         std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) = 0;

    // => Knobs <= //

    void set_early_screening(bool early_screening) { early_screening_ = early_screening; }
    void set_lr_symmetric(bool lr_symmetric) { lr_symmetric_ = lr_symmetric_; }
    /// Bench accessors
    void set_bench(int bench) { bench_ = bench; }
    int get_bench() const { return bench_; }

    /**
    * Print header information regarding JK
    * type on output file
    */
    virtual void print_header() const = 0;

    /**
    * Return number of ERI shell quartets computed during the SplitJK build process.
    */
    virtual size_t num_computed_shells();

    /**
    * print name of method
    */
    virtual std::string name() = 0;
};

// ==> Start SplitJK Coulomb (J) Algorithms here <== //

/// Build the coulomb (J) matrix using DF-DirJ
/// Reference is https://doi.org/10.1039/B204199P
class PSI_API DirectDFJ : public SplitJK {
    // => Density Fitting Stuff <= //

    /// Auxiliary basis set
    std::shared_ptr<BasisSet> auxiliary_;
    /// Coulomb Metric
    SharedMatrix J_metric_;

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
                 std::vector<std::shared_ptr<Matrix> >& G_comp,
         std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) override;

    // => Knobs <= //
    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() const override;

    /**
    * Return number of ERI shell quartets computed during the SplitJK build process.
    */
    size_t num_computed_shells() override;

    /**
    * print name of method
    */
    std::string name() override { return "DF-DirJ"; }
};

// ==> Start SplitJK Exchange (K) Algorithms here <== //

/**
 * @author Andy Jiang, Georgia Tech, December 2021
 *
 * @brief constructs the K matrix using the LinK algorithm, described in [Ochsenfeld:1998:1663]_
 * doi: 10.1063/1.476741
 */
class PSI_API LinK : public SplitJK {
    // => LinK variables <= //

    // Density-based ERI Screening tolerance to use in the LinK algorithm
    double linK_ints_cutoff_;

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    LinK(std::shared_ptr<BasisSet> primary, Options& options);
    /// Destructor
    ~LinK() override;

    /// Build the exchange (K) matrix using LinK
    void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp,
         std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) override;

    // => Knobs <= //

    /**
    * Print header information regarding SplitJK
    * type on output file
    */
    void print_header() const override;

    /**
    * Return number of ERI shell quartets computed during the SplitJK build process.
    */
    size_t num_computed_shells() override;

    /**
    * print name of method
    */
    std::string name() override { return "LinK"; }
};

/**
 * @brief constructs the K matrix using the Chain-of-Spheres Exchange (COSX)
 * algorithm, described in [Neese:2009:98]_
 * doi: 10.1016/j.chemphys.2008.10.036
 */
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

    // integral cutoff
    double kscreen_;
    // density element cutoff
    double dscreen_;
    // basis cutoff
    double basis_tol_;
    /// use overlap-fitted COSX algo?
    bool overlap_fitted_;

   public:
    // => Constructors < = //

    /**
     * @param primary primary basis set for this system.
     *        AO2USO transforms will be built with the molecule
     *        contained in this basis object, so the incoming
     *        C matrices must have the same spatial symmetry
     *        structure as this molecule
     */
    COSK(std::shared_ptr<BasisSet> primary, Options& options);
    /// Destructor
    ~COSK() override;

    /// Build the exchange (K) matrix using COSX
    /// primary reference is https://doi.org/10.1016/j.chemphys.2008.10.036
    /// overlap fitting is discussed in https://doi.org/10.1063/1.3646921
    void build_G_component(std::vector<std::shared_ptr<Matrix> >& D,
                 std::vector<std::shared_ptr<Matrix> >& G_comp,
         std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) override;

    // => Knobs <= //
    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() const override;

    /**
    * Return number of ERI shell quartets computed during the SplitJK build process.
    */
    size_t num_computed_shells() override;

    /**
    * print name of method
    */
    std::string name() override { return "COSX"; }
};

}

#endif

