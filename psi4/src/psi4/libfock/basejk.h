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

#ifndef libmints_basejk_h
#define libmints_basejk_h

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "psi4/pragma.h"

namespace psi {

class BasisSet;

/**
 * Class BaseJK
 *
 * Shared logistics base for real (JK) and complex (ComplexJK) Fock builders.
 * Holds print/debug/memory/cutoff/do_J/do_K knobs and related accessors that
 * are independent of Matrix vs ComplexMatrix storage.
 */
class PSI_API BaseJK {
   protected:
    // => Utility Variables <= //

    /// Print flag, defaults to 1
    int print_;
    /// Debug flag, defaults to 0
    int debug_;
    /// Bench flag, defaults to 0
    int bench_;
    /// Memory available, in doubles, defaults to 256 MB (32 M doubles)
    size_t memory_;
    /// Number of OpenMP threads (defaults to 1 in no OpenMP, Process::environment.get_n_threads() otherwise)
    int omp_nthread_;
    /// Integral cutoff (defaults to 0.0)
    double cutoff_;
    /// Number of ERI shell quartets computed, i.e., not screened out
    size_t num_computed_shells_;
    /// Tally of ERI shell n-lets (triplets, quartets) computed per SCF iteration
    std::unordered_map<std::string, std::vector<size_t>> computed_shells_per_iter_;

    // => Tasks <= //

    /// Do J matrices? Defaults to true
    bool do_J_;
    /// Do K matrices? Defaults to true
    bool do_K_;

    /// Left-right symmetric? Determined in each call of compute()
    bool lr_symmetric_;

    /// Primary basis set
    std::shared_ptr<BasisSet> primary_;

    /**
     * Return number of ERI shell quartets computed during the JK build process.
     */
    virtual size_t num_computed_shells();

    /// Initialize logistics knobs to defaults
    void common_init();

   public:
    /**
     * @param primary primary basis set for this system.
     */
    BaseJK(std::shared_ptr<BasisSet> primary);

    /// Destructor
    virtual ~BaseJK();

    // => Knobs <= //

    /**
     * Cutoff for individual contributions to the J/K matrices
     * Eventually we hope to use Schwarz/MBIE/Density cutoffs,
     * for now just Schwarz
     * @param cutoff ceiling of magnitude of elements to be
     *        ignored if possible
     */
    virtual void set_cutoff(double cutoff) { cutoff_ = cutoff; }
    double get_cutoff() const { return cutoff_; }
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
    void set_omp_nthread(int omp_nthread) { omp_nthread_ = omp_nthread; }
    int get_omp_nthread() const { return omp_nthread_; }

    /// Print flag (defaults to 1)
    void set_print(int print) { print_ = print; }
    /// Debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
    /// Bench flag (defaults to 0)
    void set_bench(int bench) { bench_ = bench; }
    int get_bench() const { return bench_; }
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
    virtual void set_do_K(bool do_K) { do_K_ = do_K; }

    // => Accessors <= //

    /**
     * Returns the internal primary basis set.
     */
    std::shared_ptr<BasisSet> basisset() { return primary_; }

    /**
     * Return number of ERI shell n-lets (triplets, quartets) computed per SCF iteration during the JK build process.
     */
    const std::unordered_map<std::string, std::vector<size_t>>& computed_shells_per_iter();
    const std::vector<size_t>& computed_shells_per_iter(const std::string& n_let);
};

}  // namespace psi

#endif
