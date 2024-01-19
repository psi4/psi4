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

#ifndef _psi_src_lib_libmints_twobody_h
#define _psi_src_lib_libmints_twobody_h

#include "psi4/pragma.h"

#include <functional>
#include <memory>
#include <tuple>
#include <vector>

#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/matrix.h"

namespace psi {

enum class ScreeningType { None, Schwarz, CSAM, QQR, Density };

enum PermutedOrder { ABCD = 0, BACD = 1, ABDC = 2, BADC = 3, CDAB = 4, CDBA = 5, DCAB = 6, DCBA = 7 };

typedef std::vector<std::pair<int, int>> ShellPairBlock;

class IntegralFactory;
class AOShellCombinationsIterator;
class BasisSet;
class GaussianShell;

/*! \ingroup MINTS
 *  \class TwoBodyInt
 *  \brief Two body integral base class.
 */
class PSI_API TwoBodyAOInt {
   protected:
    const IntegralFactory *integral_;

    const std::shared_ptr<BasisSet> original_bs1_;
    const std::shared_ptr<BasisSet> original_bs2_;
    const std::shared_ptr<BasisSet> original_bs3_;
    const std::shared_ptr<BasisSet> original_bs4_;

    std::shared_ptr<BasisSet> bs1_;
    std::shared_ptr<BasisSet> bs2_;
    std::shared_ptr<BasisSet> bs3_;
    std::shared_ptr<BasisSet> bs4_;

    /// Buffer to hold the final integrals.
    double *target_full_;

    /// Pointers to each chunk of derivative integrals
    std::vector<const double *> buffers_;

    /// Where to put the next integrals (should be part of target_full_)
    double *target_;

    /// Number of integrals in the current buffer
    int curr_buff_size_;
    /// Buffer to hold the transformation intermediates.
    double *tformbuf_;
    /// Buffer to hold the initially computed integrals.
    double *source_full_;
    /// Where to put the next temporary integrals (should be part of source_full_)
    double *source_;
    /// Maximum number of unique quartets needed to compute a set of SO's
    int max_unique_quartets_;
    /// Number of atoms.
    int natom_;
    /// Derivative level.
    int deriv_;
    /// How the shells were reordered for libint
    PermutedOrder permuted_order_;
    /// Are the basis sets in the bra the same?
    bool bra_same_;
    /// Are the basis sets in the ket the same?
    bool ket_same_;
    /// Are the basis sets in the bra and the ket all the same?
    bool braket_same_;

    /// The blocking scheme used for the integrals
    std::vector<ShellPairBlock> blocks12_, blocks34_;

    /*
     * Sieve information
     */
    typedef std::vector<std::pair<int, int>> PairList;
    /// The threshold below which integrals are to be neglected
    double screening_threshold_;
    double screening_threshold_squared_;
    int nshell_;
    int nbf_;
    /// The algorithm to use for screening
    ScreeningType screening_type_;
    /// |(mn|mn)| values (nbf * nbf)
    std::vector<double> function_pair_values_;
    /// max |(MN|MN)| values (nshell * nshell)
    std::vector<double> shell_pair_values_;
    /// max |(MM|NN)| values (nshell * nshell)
    std::vector<double> shell_pair_exchange_values_;
    /// sqrt|(mm|mm)| values (nshell)
    std::vector<double> function_sqrt_;
    /// Max density per matrix (Outer loop over density matrices, inner loop over shell pairs)
    std::vector<std::vector<double>> max_dens_shell_pair_;
    /// Significant unique function pairs, in row-major, lower triangular indexing
    PairList function_pairs_;
    /// Significant unique shell pairs, in row-major, lower triangular indexing
    PairList shell_pairs_, shell_pairs_bra_, shell_pairs_ket_;
    /// The largest value of any integral as predicted by the sieving method
    double max_integral_;
    /// Unique function pair indexing, accessed in triangular order, or -1 for non-significant pair
    std::vector<long int> function_pairs_reverse_;
    /// Unique shell pair indexing, accessed in triangular order, or -1 for non-significant pair
    std::vector<long int> shell_pairs_reverse_;
    /// Significant function pairs, indexes by function
    std::vector<std::vector<int>> shell_to_shell_;
    /// Significant shell pairs, indexes by shell
    std::vector<std::vector<int>> function_to_function_;
    std::function<bool(int, int, int, int)> sieve_impl_;

    void setup_sieve();
    void create_sieve_pair_info(const std::shared_ptr<BasisSet> bs, PairList &shell_pairs, bool is_bra);

    /// Implements CSAM screening of a shell quartet
    bool shell_significant_csam(int M, int N, int R, int S);
    /// Implements Schwarz inequality screening of a shell quartet
    bool shell_significant_schwarz(int M, int N, int R, int S);
    /// Asks whether this shell quartet contributes by the density test (Haser 1989)
    bool shell_significant_density(int M, int N, int R, int S);
    /// Implements the null screening of a shell quartet - always true
    bool shell_significant_none(int M, int N, int R, int S);

    /*! Create the optimal blocks of shell pairs
     *
     * Default implementation
     */
    virtual void create_blocks();

    void permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34, bool p13p24);
    void permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);
    void permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4);

    //    TwoBodyInt(std::shared_ptr<BasisSet> bs1,
    //               std::shared_ptr<BasisSet> bs2,
    //               std::shared_ptr<BasisSet> bs3,
    //               std::shared_ptr<BasisSet> bs4,
    //               int deriv = 0);

    TwoBodyAOInt(const IntegralFactory *intsfactory, int deriv = 0);

    TwoBodyAOInt(const TwoBodyAOInt &rhs);

    virtual size_t compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int sh1, int sh2, int sh3, int sh4, bool is_bra) = 0;
   public:
    virtual ~TwoBodyAOInt();

    /// Basis set on center one
    std::shared_ptr<BasisSet> basis();
    /// Basis set on center one
    std::shared_ptr<BasisSet> basis1();
    /// Basis set on center two
    std::shared_ptr<BasisSet> basis2();
    /// Basis set on center three
    std::shared_ptr<BasisSet> basis3();
    /// Basis set on center four
    std::shared_ptr<BasisSet> basis4();

    /*
     * Sieve information
     */
    /// Update max_dens_shell_pair_ given an updated density matrix (Haser 1989)
    void update_density(const std::vector<SharedMatrix>& D);
    /// Ask the built in sieve whether this quartet contributes
    bool shell_significant(int M, int N, int R, int S) const { return sieve_impl_(M, N, R, S); };
    /// Are any of the quartets within a given shellpair list significant
    bool shell_block_significant(int shellpair12, int shellpair34) const;
    /// Does a given shell pair contribute to any significant integrals?
    bool shell_pair_significant(int shell1, int shell2) const;
    /// Square of ceiling of shell quartet (MN|RS)
     inline double shell_ceiling2(int M, int N, int R, int S) {
        return shell_pair_values_[N * nshell_ + M] * shell_pair_values_[R * nshell_ + S];
    }
    /// Is the function pair (mn| ever significant according to sieve (no restriction on mn order)
    inline bool function_pair_significant(const int m, const int n) {
        return function_pair_values_[m * nbf_ + n] * max_integral_ >= screening_threshold_squared_;
    }
    /// Is the integral (mn|rs) significant according to sieve? (no restriction on mnrs order)
    inline bool function_significant(const int m, const int n, const int r, const int s) {
        return function_pair_values_[m * nbf_ + n] * function_pair_values_[r * nbf_ + s] >= screening_threshold_squared_;
    }
    /// Return max(PQ|PQ)
    double max_integral() const { return max_integral_; }
    /// Square of ceiling of integral (mn|rs)
     inline double function_ceiling2(int m, int n, int r, int s) {
        return function_pair_values_[m * nbf_ + n] * function_pair_values_[r * nbf_ + s];
    }
    // the value of the bound for pair m and n
    double shell_pair_value(int m, int n) { return shell_pair_values_[m * nshell_ + n]; };
    /// Return the maximum density matrix element per shell pair. Maximum is over density matrices, if multiple set
    double shell_pair_max_density(int M, int N) const;

    /// For a given PQ shellpair index, what's the first RS pair that should be processed such
    /// that loops may be processed generating only permutationally unique PQ<=RS.  For engines
    /// that don't use blocking in the ket it is trivial, but a little more complex with ket blocks
    virtual size_t first_RS_shell_block(size_t PQpair) const { return PQpair; }

    /// Significant unique function pair list, with only m>=n elements listed
    const std::vector<std::pair<int, int> >& function_pairs() const { return function_pairs_; }
    /// Significant unique shell pair pair list, with only M>=N elements listed
    const std::vector<std::pair<int, int> >& shell_pairs() const { return shell_pairs_; }
    /// Unique function pair indexing, element m*(m+1)/2 + n (where m>=n) gives the dense index or
    /// -1 if the function pair does not contribute
    const std::vector<long int> function_pairs_to_dense() const { return function_pairs_reverse_; }
    /// Unique shell pair indexing, element M*(M+1)/2 + N (where M>=N) gives the dense index or
    /// -1 if the shell pair does not contribute
    const std::vector<long int> shell_pairs_to_dense() const { return shell_pairs_reverse_; }
    /// Significant function pairs; for each function it gives a list of functions that contribute to make a function pair
    const std::vector<std::vector<int> >& significant_partners_per_function() const { return function_to_function_; }
    /// Significant shell pairs; for each shell it gives a list of shells that contribute to make a shell pair
    const std::vector<std::vector<int> >& significant_parterns_per_shell() const { return shell_to_shell_; }

    /// Returns the derivative level this object is setup for.
    int deriv() const { return deriv_; }

    /// Buffer where the integrals are placed
    const double *buffer() const { return target_full_; }

    /// Buffer where each chunk of integrals is placed
    const std::vector<const double *> &buffers() const { return buffers_; }

    /// The maximum number of batched functions in the ket (1 for a non vectorized engine).
    virtual int maximum_block_size() const { return 1; }

    /// Returns the integral factory used to create this object
    const IntegralFactory *integral() const { return integral_; }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual size_t compute_shell(const AOShellCombinationsIterator &) = 0;

    //! Get optimal blocks of shell pairs for centers 1 & 2
    std::vector<ShellPairBlock> get_blocks12() const;

    //! Get optimal blocks of shell pairs for centers 3 & 4
    std::vector<ShellPairBlock> get_blocks34() const;

    /*! Compute integrals for two blocks
     *
     * The indices \p shellpair12 and \p shellpair34 refer to the indices
     * in the vectors returned by get_blocks12() and get_blocks34(),
     * respectively.
     *
     * The parameters \p npair12 and \p npair34 refer to how many
     * of the shell pairs of centers 1 & 2 and 3 & 4, respectively,
     * to actually calculate. This can be used to compute only
     * triangular parts, for example.  A value of -1 means to calculate
     * all that are part of the shell pair batch.
     */
    virtual void compute_shell_blocks(int shellpair12, int shellpair34, int npair12 = -1, int npair34 = -1);
    /*! Compute derivative integrals for two blocks */
    virtual void compute_shell_blocks_deriv1(int shellpair12, int shellpair34, int npair12 = -1, int npair34 = -1);
    /*! Compute derivative integrals for two blocks */
    virtual void compute_shell_blocks_deriv2(int shellpair12, int shellpair34, int npair12 = -1, int npair34 = -1);

    /// Is the shell zero?
    virtual int shell_is_zero(int, int, int, int) { return 0; }

    virtual size_t compute_shell(int s1, int s2, int s3, int s4) = 0;

    /// Compute the first derivatives
    virtual size_t compute_shell_deriv1(int s1, int s2, int s3, int s4) = 0;

    /// Compute the second derivatives
    virtual size_t compute_shell_deriv2(int s1, int s2, int s3, int s4) = 0;

    /// Normalize Cartesian functions based on angular momentum
    void normalize_am(std::shared_ptr<GaussianShell>, std::shared_ptr<GaussianShell>, std::shared_ptr<GaussianShell>,
                      std::shared_ptr<GaussianShell>, int nchunk = 1);

    /// Returns a clone of this object. By default throws an exception
    virtual TwoBodyAOInt *clone()  const = 0;

    /// Results go back to buffer_
    void pure_transform(int, int, int, int, int nchunk, bool copy_to_source = true);
};

typedef std::shared_ptr<TwoBodyAOInt> SharedTwoBodyAOInt;

}  // namespace psi

#endif
