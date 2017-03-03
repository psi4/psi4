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

#ifndef JK_GRAD_H
#define JK_GRAD_H

#include "psi4/libmints/typedefs.h"
#include <map>

namespace psi {

class ERISieve;

namespace scfgrad {

class JKGrad {

protected:
    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Memory available, in doubles, defaults to 256 MB (32 M doubles)
    unsigned long int memory_;
    /// Number of OpenMP threads (defaults to 1 in no OpenMP, Process::environment.get_n_threads() otherwise)
    int omp_num_threads_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;
    /// Maximum derivative level
    int deriv_;

    std::shared_ptr<BasisSet> primary_;

    /// Sieve, must be static throughout the life of the object
    std::shared_ptr<ERISieve> sieve_;

    SharedMatrix Ca_;
    SharedMatrix Cb_;

    SharedMatrix Da_;
    SharedMatrix Db_;
    SharedMatrix Dt_;

    bool do_J_;
    bool do_K_;
    bool do_wK_;

    double omega_;

    std::map<std::string, SharedMatrix> gradients_;
    std::map<std::string, SharedMatrix> hessians_;

    void common_init();

public:
    JKGrad(int deriv, std::shared_ptr<BasisSet> primary);
    virtual ~JKGrad();

    /**
    * Static instance constructor, used to get prebuilt DFJK/DirectJK objects
    * using knobs in options.
    * @param options Options reference, with preset parameters
    * @return abstract JK object, tuned in with preset options
    */
    static std::shared_ptr<JKGrad> build_JKGrad(int deriv, std::shared_ptr<BasisSet> primary,
                                                std::shared_ptr<BasisSet> auxiliary);

    void set_Ca(SharedMatrix Ca) { Ca_ = Ca; }
    void set_Cb(SharedMatrix Cb) { Cb_ = Cb; }
    void set_Da(SharedMatrix Da) { Da_ = Da; }
    void set_Db(SharedMatrix Db) { Db_ = Db; }
    void set_Dt(SharedMatrix Dt) { Dt_ = Dt; }

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
    void set_memory(unsigned long int memory) { memory_ = memory; }
    /**
     * Maximum number of OpenMP threads to use. It may be necessary
     * to clamp this to some value smaller than the total number of
     * cores for machines with a high core-to-memory ratio to avoid
     * running out of memory due to integral generation objects
     * @param omp_nthread Maximum number of threads to use in
     *        integral generation objects (BLAS/LAPACK can still
     *        run with their original maximum number)
     */
    void set_omp_num_threads(int omp_nthread) { omp_num_threads_ = omp_nthread; }
    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }
    /**
    * Set to do J tasks
    * @param do_J do J matrices or not,
    *        defaults to true
    */
    void set_do_J(bool do_J) { do_J_ = do_J; }
    /**
    * Set to do K tasks
    * @param do_K do K matrices or not,
    *        defaults to true
    */
    void set_do_K(bool do_K) { do_K_ = do_K; }
    /**
    * Set to do wK tasks
    * @param do_wK do wK matrices or not,
    *        defaults to false
    */
    void set_do_wK(bool do_wK) { do_wK_ = do_wK; }
    /**
    * Set the omega value for wK
    * @param omega range-separation parameter
    */
    void set_omega(double omega) { omega_ = omega; }

    std::map<std::string, SharedMatrix>& gradients() { return gradients_; }
    std::map<std::string, SharedMatrix>& hessians() { return hessians_; }

    virtual void compute_gradient() = 0;
    virtual void compute_hessian() = 0;

    virtual void print_header() const = 0;
};

class DFJKGrad : public JKGrad {

protected:
    std::shared_ptr<BasisSet> auxiliary_;

    std::shared_ptr<PSIO> psio_;

    /// Number of threads for DF integrals
    int df_ints_num_threads_;
    /// Condition cutoff in fitting metric, defaults to 1.0E-12
    double condition_;

    void common_init();

    void build_Amn_terms();
    void build_Amn_lr_terms();
    void build_AB_inv_terms();
    void build_UV_terms();
    void build_AB_x_terms();
    void build_Amn_x_terms();
    void build_Amn_x_lr_terms();

    /// File number for Alpha (Q|mn) tensor
    unsigned int unit_a_;
    /// File number for Beta (Q|mn) tensor
    unsigned int unit_b_;
    /// File number for J tensors
    unsigned int unit_c_;

public:
    DFJKGrad(int deriv, std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary);
    virtual ~DFJKGrad();

    void compute_gradient();
    void compute_hessian();

    void print_header() const;

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
    void set_unit_a(unsigned int unit) { unit_a_ = unit; }
    /**
     * Which file number should the Beta (Q|mn) integrals go in
     * @param unit Unit number
     */
    void set_unit_b(unsigned int unit) { unit_b_ = unit; }
    /**
     * Which file number should the J tensors go in
     * @param unit Unit number
     */
    void set_unit_c(unsigned int unit) { unit_c_ = unit; }

    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_df_ints_num_threads(int val) { df_ints_num_threads_ = val; }
};

class DirectJKGrad : public JKGrad {

protected:
    // Number of threads to use
    int ints_num_threads_;

    void common_init();

    std::map<std::string, std::shared_ptr<Matrix> > compute1(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints);
    std::map<std::string, std::shared_ptr<Matrix> > compute2(std::vector<std::shared_ptr<TwoBodyAOInt> >& ints);
public:
    DirectJKGrad(int deriv, std::shared_ptr<BasisSet> primary);
    virtual ~DirectJKGrad();

    void compute_gradient();
    void compute_hessian();

    void print_header() const;

    /**
     * What number of threads to compute integrals on
     * @param val a positive integer
     */
    void set_ints_num_threads(int val) { ints_num_threads_ = val; }


};

}} // Namespaces
#endif
