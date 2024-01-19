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

#ifndef SIEVE_H
#define SIEVE_H

// need this for erfc^{-1} in the QQR sieve
//#include <cfloat>
#include <vector>
#include <memory>
//#include <utility>
#include "psi4/pragma.h"
#include "psi4/libmints/vector3.h"

namespace psi {

class BasisSet;

/**
 * ERISieve
 *
 * Class to perform ERI sieving
 *
 * At the moment, this class uses Cauchy-Schwarz sieving,
 * and is designed to be used in concert with density-based
 * sieving. In the future, it is hoped that MBIE sieving may
 * be added to this class.
 *
 * Use of this class is as follows:
 *
 *     // Initialize the sieve object
 *     std::shared_ptr<ERISieve> sieve(basisset, sieve_cutoff);
 *
 *     // Reset the sieve cutoff (you can do this wherever)
 *     sieve->set_sieve(new_cutoff);
 *
 *     ...
 *
 *     // Investigate stuff. This example is for shells, all methods are
 *     // also provided for functions
 *
 *     // Compute a shell quartet (MN|RS), if (MN|RS) >= sieve_cutoff
 *     if (sieve->shell_significant(M,N,R,S)) eri->compute(M,N,R,S);
 *
 *     // Compute a shell quartet (MN|RS), if (MN|RS) * D_RS >= sieve_cutoff
 *     // Squares are used to avoid sqrt()
 *     if (sieve->shell_ceiling2(M,N,R,S) * D_RS * D_RS >= sieve_cutoff * sieve_cutoff)
 *         eri->compute(M,N,R,S);
 *
 *     // Index the significant MN shell pairs (triangular M,N)
 *     const std::vector<std::pair<int,int> >& MN = sieve->shell_pairs();
 *     for (long int index = 0L; index < MN.size(); ++index) {
 *         int M = MN[index].first;
 *         int N = MN[index].second;
 *     }
 *
 *     // Check if a triangular index MNindex (M * (M + 1) / 2) + N exists,
 *     // and if so, where it starts in reduced triangular MN
 *     int MNindex = (M * (M + 1) >> 1) + N;
 *     const std::vector<long int> >& MN_reverse = sieve->shell_pairs_reverse();
 *     int MNreduced = MN_reverse[MNindex];
 *     if (MNreduced < 0) {
 *         // The shell pair is not signficant
 *     } else {
 *         // The shell pair is the MNreduced significant shell pair
 *     }
 *
 *
 */

class PSI_API __attribute__((deprecated(
    "ERISieve is deprecated in favor of TwoBodyAOInt, and "
    "will be fully be removed as soon as Psi4 v1.9 releases. "
))) ERISieve {
   protected:
    /// Debug flag (defaults to 0)
    int debug_;

    /// Basis set reference
    std::shared_ptr<BasisSet> primary_;

    /// Number of basis functions
    size_t nbf_;
    /// Number of shells
    size_t nshell_;

    /// Cutoff values
    double sieve_;
    /// Maximum |(mn|ls)|
    double max_;
    /// sieve_ / max_
    double sieve_over_max_;
    /// sieve_ * sieve_
    double sieve2_;
    /// sieve_ * sieve_ / max_
    double sieve2_over_max_;

    /// |(mn|mn)| values (nbf * nbf)
    std::vector<double> function_pair_values_;
    /// max |(MN|MN)| values (nshell * nshell)
    std::vector<double> shell_pair_values_;

    /// Significant unique bra- function pairs, in reduced triangular indexing
    std::vector<std::pair<int, int> > function_pairs_;
    /// Significant unique bra- shell pairs, in reduced triangular indexing
    std::vector<std::pair<int, int> > shell_pairs_;
    /// Unique bra- function pair indexing, accessed in triangular order, or -1 for non-significant pair
    std::vector<long int> function_pairs_reverse_;
    /// Unique bra- shell pair indexing, accessed in triangular order, or -1 for non-significant pair
    std::vector<long int> shell_pairs_reverse_;
    /// Significant function pairs, indexes by function
    std::vector<std::vector<int> > shell_to_shell_;
    /// Significant shell pairs, indexes by shell
    std::vector<std::vector<int> > function_to_function_;

    ///////////////////////////////////////
    // adding stuff for QQR sieves

    bool do_qqr_;

    // erfc^{-1}(threshold), used in QQR sieving
    double erfc_thresh_;

    // need an array of extents from the definition

    // key: how do I efficiently check integrals? without conditional on which
    // screening I do
    //
    // 1) just a different function called outside - shell_significant_qqr()

    // integrals() - fills array of extents

    // r_{\mu \nu} in QQR paper (eqn. B2)
    std::vector<Vector3> contracted_centers_;

    // ext'_{\mu \nu} (Eqn. B4)
    // Extents of contracted charge distributions
    std::vector<double> extents_;

    ////////////////////////////////////////
    // adding stuff for CSAM sieving (DOI 10.1063/1.4994190)

    bool do_csam_;

    /// max |(MM|NN)| values (nshell * nshell)
    std::vector<double> shell_pair_exchange_values_;
    /// sqrt|(mm|mm)| values (nshell)
    std::vector<double> function_sqrt_;
    /// Compute csam sieve integrals (only done once)
    void csam_integrals();

    ///////////////////////////////////////

    /// Set initial indexing
    void common_init();
    /// Compute sieve integrals (only done once)
    void integrals();

   public:
    /// Constructor, basis set and first sieve cutoff
    ERISieve(std::shared_ptr<BasisSet> primary, double sieve = 0.0, bool do_csam = false);
    /// Destructor, frees memory
    virtual ~ERISieve();

    /// Set sieve value and redo indexing
    void set_sieve(double sieve);
    /// Get sieve cutoff value
    double sieve() const { return sieve_; }
    /// Global maximum |(mn|rs)|
    double max() const { return max_; }
    /// Get whether this sieve performs CSAM screening (as opposed to vanilla Schwarz)
    bool do_csam() const { return do_csam_; }

    // => Significance Checks <= //

    /// Square of ceiling of shell quartet (MN|RS)
    inline double shell_ceiling2(int M, int N, int R, int S) {
        return shell_pair_values_[N * nshell_ + M] * shell_pair_values_[R * nshell_ + S];
    }

    /// Square of ceiling of integral (mn|rs)
    inline double function_ceiling2(int m, int n, int r, int s) {
        return function_pair_values_[m * nbf_ + n] * function_pair_values_[r * nbf_ + s];
    }

    /// Is the shell quartet (MN|RS) significant according to sieve? (no restriction on MNRS order)

    // inline bool shell_significant(int M, int N, int R, int S) {
    bool shell_significant(int M, int N, int R, int S) {
        bool schwarz_bound =
            shell_pair_values_[N * nshell_ + M] * shell_pair_values_[R * nshell_ + S] >= sieve2_;
        if (do_qqr_ && schwarz_bound) {
            return shell_significant_qqr(M, N, R, S);
        } else if (do_csam_ && schwarz_bound) {
            return shell_significant_csam(M, N, R, S);
        } else {
            return schwarz_bound;
        }
    }

    // Implements the QQR sieve
    bool shell_significant_qqr(int M, int N, int R, int S);

    // Implements the CSAM sieve
    bool shell_significant_csam(int M, int N, int R, int S);

    /// Is the integral (mn|rs) significant according to sieve? (no restriction on mnrs order)
    inline bool function_significant(int m, int n, int r, int s) {
        return function_pair_values_[m * nbf_ + n] * function_pair_values_[r * nbf_ + s] >= sieve2_;
    }

    /// Is the shell pair (MN| ever significant according to sieve (no restriction on MN order)
    inline bool shell_pair_significant(int M, int N) {
        return shell_pair_values_[M * nshell_ + N] * max_ >= sieve2_;
    }

    /// Is the function pair (mn| ever significant according to sieve (no restriction on mn order)
    inline bool function_pair_significant(int m, int n) {
        return function_pair_values_[m * nbf_ + n] * max_ >= sieve2_;
    }
    // => Indexing [these change after a call to sieve()] <= //


    // void shell_pair_values(std::vector<std::vector<std::pair<double, int> > >& values) const;

    // just return the value of the bound for pair m and n
    double shell_pair_value(int m, int n) const;
    // return the vector of
    std::vector<double> shell_pair_values() { return shell_pair_values_; }
    // return the vector of function pairs
    std::vector<double> function_pair_values() { return function_pair_values_; }

    /// Significant unique bra- function pairs, in reduced triangular indexing
    const std::vector<std::pair<int, int> >& function_pairs() const { return function_pairs_; }
    /// Significant unique bra- shell pairs, in reduced triangular indexing
    const std::vector<std::pair<int, int> >& shell_pairs() const { return shell_pairs_; }
    /// Unique bra- function pair indexing, accessed in triangular order, or -1 for non-significant pair
    const std::vector<long int> function_pairs_reverse() const { return function_pairs_reverse_; }
    /// Unique bra- shell pair indexing, accessed in triangular order, or -1 for non-significant pair
    const std::vector<long int> shell_pairs_reverse() const { return shell_pairs_reverse_; }
    /// Significant function pairs, indexes by function
    const std::vector<std::vector<int> >& function_to_function() const { return function_to_function_; }
    /// Significant shell pairs, indexes by shell
    const std::vector<std::vector<int> >& shell_to_shell() const { return shell_to_shell_; }

    /// Set debug flag (defaults to 0)
    void set_debug(int debug) { debug_ = debug; }
};

}  // namespace psi
#endif
