/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#ifndef Corr_GRAD_H
#define Corr_GRAD_H

#include "psi4/libmints/typedefs.h"
#include <map>

namespace psi {

class ERISieve;
class BasisSet;
class PSIO;

namespace dfmp2 {

class CorrGrad {
   protected:
    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Memory available, in doubles, defaults to 256 MB (32 M doubles)
    size_t memory_;
    /// Number of OpenMP threads (defaults to 1 in no OpenMP, Process::environment.get_n_threads() otherwise)
    int nthreads_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;

    std::shared_ptr<BasisSet> primary_;

    /**
     * Rules:
     * CC' = D (spin)
     * LL' - RR' = P (spin)
     * Da + Db = D
     * Pa + Pb = P
     **/

    SharedMatrix Ca_;
    SharedMatrix Cb_;
    SharedMatrix La_;
    SharedMatrix Lb_;
    SharedMatrix Ra_;
    SharedMatrix Rb_;

    SharedMatrix Da_;
    SharedMatrix Db_;
    SharedMatrix Dt_;
    SharedMatrix Pa_;
    SharedMatrix Pb_;
    SharedMatrix Pt_;

    std::map<std::string, SharedMatrix> gradients_;

    void common_init();

   public:
    CorrGrad(std::shared_ptr<BasisSet> primary);
    virtual ~CorrGrad();

    /**
     * Static instance constructor, used to get prebuilt DFCorr/DirectCorr objects
     * using knobs in options.
     * @param options Options reference, with preset parameters
     * @return abstract Corr object, tuned in with preset options
     */
    static std::shared_ptr<CorrGrad> build_CorrGrad(std::shared_ptr<BasisSet> primary,
                                                    std::shared_ptr<BasisSet> auxiliary);

    void set_Ca(SharedMatrix Ca) { Ca_ = Ca; }
    void set_Cb(SharedMatrix Cb) { Cb_ = Cb; }
    void set_La(SharedMatrix La) { La_ = La; }
    void set_Lb(SharedMatrix Lb) { Lb_ = Lb; }
    void set_Ra(SharedMatrix Ra) { Ra_ = Ra; }
    void set_Rb(SharedMatrix Rb) { Rb_ = Rb; }
    void set_Da(SharedMatrix Da) { Da_ = Da; }
    void set_Db(SharedMatrix Db) { Db_ = Db; }
    void set_Dt(SharedMatrix Dt) { Dt_ = Dt; }
    void set_Pa(SharedMatrix Pa) { Pa_ = Pa; }
    void set_Pb(SharedMatrix Pb) { Pb_ = Pb; }
    void set_Pt(SharedMatrix Pt) { Pt_ = Pt; }

    /**
     * Cutoff for individual contributions to the J/K matrices
     * Eventually we hope to use Schwarz/MBIE/Density cutoffs,
     * for now just Schwarz
     * @param cutoff ceiling of magnitude of elements to be
     *        ignored if possible
     */
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }
    /**
     * Maximum memory to use, in doubles (for tensor-based methods,
     * integral generation objects typically ignore this)
     * @param memory maximum number of doubles to allocate
     */
    void set_memory(size_t memory) { memory_ = memory; }
    /**
     * Maximum number of OpenMP threads to use. It may be necessary
     * to clamp this to some value smaller than the total number of
     * cores for machines with a high core-to-memory ratio to avoid
     * running out of memory due to integral generation objects
     * @param omp_nthread Maximum number of threads to use in
     *        integral generation objects (BLAS/LAPACK can still
     *        run with their original maximum number)
     */
    void set_omp_num_threads(int omp_nthread) { nthreads_ = omp_nthread; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }

    std::map<std::string, SharedMatrix>& gradients() { return gradients_; }

    virtual void compute_gradient() = 0;

    virtual void print_header() const = 0;
};

class DFCorrGrad : public CorrGrad {
   protected:
    std::shared_ptr<BasisSet> auxiliary_;

    std::shared_ptr<PSIO> psio_;

    /// Number of threads for DF integrals
    int df_ints_num_threads_;
    /// Condition cutoff in fitting metric, defaults to 1.0E-12
    double condition_;

    /// Sieve, must be static throughout the life of the object
    std::shared_ptr<ERISieve> sieve_;

    void common_init();

    void build_Amn_terms();
    void build_AB_inv_terms();
    void build_UV_terms();
    void build_AB_x_terms();
    void build_Amn_x_terms();

    void fitting_helper(SharedMatrix J, size_t file, const std::string& label, size_t naux, size_t nij, size_t memory);
    void UV_helper(SharedMatrix V, double c, size_t file, const std::string& label, size_t naux, size_t nij,
                   size_t memory);

    /// File number for Alpha (Q|mn) tensor
    size_t unit_a_;
    /// File number for Beta (Q|mn) tensor
    size_t unit_b_;
    /// File number for J tensors
    size_t unit_c_;

   public:
    DFCorrGrad(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary);
    ~DFCorrGrad() override;

    void compute_gradient() override;

    void print_header() const override;

    /**
     * Minimum relative eigenvalue to retain in fitting inverse
     * All eigenvectors with \epsilon_i < condition * \epsilon_max
     * will be discarded
     * @param condition, minimum relative eigenvalue allowed,
     *        defaults to 1.0E-12
     */
    void set_condition(double condition) { condition_ = condition; }
    /**
     * Which file number should the Alpha (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit_a(size_t unit) { unit_a_ = unit; }
    /**
     * Which file number should the Beta (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit_b(size_t unit) { unit_b_ = unit; }
    /**
     * Which file number should the J tensors go in
     * @param unit Unit number
     */
    void set_unit_c(size_t unit) { unit_c_ = unit; }

    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }
};

}  // namespace dfmp2
}  // namespace psi
#endif
