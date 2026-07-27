/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef COMPLEXJK_H
#define COMPLEXJK_H

#include <memory>
#include <string>
#include <vector>

#include "psi4/pragma.h"
#include "psi4/libfock/basejk.h"
#include "psi4/libmints/complexwavefunction.h"
#include "psi4/libmints/typedefs.h"

namespace psi {

class BasisSet;
class Options;

/**
 * Class ComplexJK
 *
 * Complex-orbital counterpart to JK for Generalized Fock-Matrix Contributions
 * of the form:
 *
 *  J_mn = (mn|ls) C_li^left C_si^right
 *  K_mn = (ml|ns) C_li^left C_si^right
 *
 * API intentionally mirrors JK (push C_left, compute(), grab D/J/K), with
 * SharedComplexMatrix storage instead of SharedMatrix.
 *
 * Intentional omissions:
 *   - DFT/range-separated quantities (wK, omega, ...) are intentionally omitted.
 *   - Symmetry: C1 is explicitly enforced.
 *      - `*_ao_` varibles.
 *
 * Method bodies that require algorithm or tensor logic currently throw
 * pybind11::attribute_error until implemented.
 */
class PSI_API ComplexJK : public BaseJK {
   protected:
    // => Utility Variables <= //

    /// CSAM Screening (defaults to false)
    double do_csam_;
    /// Whether to allow desymmetrization, for cases when it's already been performed elsewhere
    std::vector<bool> input_symmetry_cast_map_;

    // => Architecture-Level State Variables (Spatial Symmetry) <= //

    /// Pseudo-occupied C matrices, left side
    std::vector<SharedComplexMatrix> C_left_;
    /// Pseudo-occupied C matrices, right side
    std::vector<SharedComplexMatrix> C_right_;
    /// Pseudo-density matrices \f$D_{ls}=C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedComplexMatrix> D_;
    /// J matrices: \f$J_{mn}=(mn|ls)C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedComplexMatrix> J_;
    /// K matrices: \f$K_{mn}=(ml|ns)C_{li}^{left}C_{si}^{right}\f$
    std::vector<SharedComplexMatrix> K_;

    // => Per-Iteration Setup/Finalize Routines <= //

    /// Build the pseudo-density D_, before compute_JK()
    void compute_D() override;
    /// Allocate J_/K_
    void allocate_JK() override;

    /// Function that sets a number of flags and allocates memory.
    void common_init() override;

    // => Helper Routines <= //

    /// Memory (doubles) used to hold J/K/C/D and ao versions, at current moment
    size_t memory_overhead() const;
    /// Zero out all J and K matrices
    void zero();

   public:
    // => Constructors <= //

    /**
     * @param primary primary basis set for this system.
     */
    ComplexJK(std::shared_ptr<BasisSet> primary);

    /// Destructor
    ~ComplexJK() override;

    /**
     * Static instance constructor, used to get prebuilt ComplexJK objects
     * using knobs in options.
     */
    static std::shared_ptr<ComplexJK> build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                               Options& options);
    static std::shared_ptr<ComplexJK> build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                               Options& options, std::string jk_type);
    static std::shared_ptr<ComplexJK> build_JK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary,
                                               Options& options, bool do_wK, size_t doubles);

    // => Knobs <= //

    /**
     * @param do_csam whether to perform CSAM screening instead of
     *      classic Schwarz screening
     */
    void set_csam(bool do_csam) { do_csam_ = do_csam; }
    double get_csam() const { return do_csam_; }

    // => Computers <= //

    /**
     * Compute D/J/K for the current C
     * Update values in your reference to
     * C_left/C_right BEFORE calling this,
     * renew your references to the matrices
     * in D/J/K AFTER calling this.
     */
    void compute() override;

    // => Accessors <= //

    /**
     * Reference to C_left queue. It is YOUR job to
     * allocate and fill this object out
     */
    std::vector<SharedComplexMatrix>& C_left() { return C_left_; }
    /**
     * Reference to C_right queue. It is YOUR job to
     * allocate and fill this object out. Only fill
     * C_left if symmetric.
     */
    std::vector<SharedComplexMatrix>& C_right() { return C_right_; }

    /**
     * Reference to J results. The reference to the
     * std::vector<SharedComplexMatrix> is valid
     * throughout the life of the object. However, the
     * entries may be changed in each call of compute();
     * @return J vector of J matrices
     */
    const std::vector<SharedComplexMatrix>& J() const { return J_; }
    /**
     * Reference to K results. The reference to the
     * std::vector<SharedComplexMatrix> is valid
     * throughout the life of the object. However, the
     * entries may be changed in each call of compute();
     * @return K vector of K matrices
     */
    const std::vector<SharedComplexMatrix>& K() const { return K_; }
    /**
     * Reference to D results. The reference to the
     * std::vector<SharedComplexMatrix> is valid
     * throughout the life of the object. However, the
     * entries may be changed in each call of compute();
     * @return D vector of D matrices
     */
    const std::vector<SharedComplexMatrix>& D() const { return D_; }
};

class PSI_API ComplexDirectJK : public ComplexJK {
    std::string name() override { return "ComplexDirectJK"; }
    size_t memory_estimate() override;

    /// Options object
    Options& options_;

    /// Do we need to backtransform to C1 under the hood?
    bool C1() const override { return false; }
    /// Setup integrals, files, etc
    void preiterations() override;
    /// Compute J/K for current C/D
    void compute_JK() override;
    /// Delete integrals, files, etc
    void postiterations() override;

    /// Common initialization
    void common_init() override;

  public:
    // => Constructors < = //
    ComplexDirectJK(std::shared_ptr<BasisSet> primary, Options& options);
    /// Destructor
    ~ComplexDirectJK() override;

    // => Accessors <= //

    /**
    * Print header information regarding JK
    * type on output file
    */
    void print_header() const override;
};

}  // namespace psi

#endif
