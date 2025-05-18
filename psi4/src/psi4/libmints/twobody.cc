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

#include <algorithm>
#include <stdexcept>
#include "psi4/libqt/qt.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libpsi4util/process.h"

#include "libint2/shell.h"

using namespace psi;

static void transform2e_1(int, SphericalTransformIter &, double *, double *, int);
static void transform2e_2(int, SphericalTransformIter &, double *, double *, int, int, int);
static void transform2e_3(int, SphericalTransformIter &, double *, double *, int, int, int);
static void transform2e_4(int, SphericalTransformIter &, double *, double *, int, int);

TwoBodyAOInt::TwoBodyAOInt(const IntegralFactory *intsfactory, int deriv)
    : integral_(intsfactory),
      original_bs1_(integral_->basis1()),
      original_bs2_(integral_->basis2()),
      original_bs3_(integral_->basis3()),
      original_bs4_(integral_->basis4()),
      bs1_(original_bs1_),
      bs2_(original_bs2_),
      bs3_(original_bs3_),
      bs4_(original_bs4_),
      target_full_(nullptr),
      target_(nullptr),
      source_full_(nullptr),
      source_(nullptr),
      deriv_(deriv) {
    // The derived classes allocate this memory.
    tformbuf_ = nullptr;
    natom_ = original_bs1_->molecule()->natom();  // This assumes the 4 bases come from the same molecule.

    // Figure out how equivalent
    bra_same_ = (original_bs1_ == original_bs2_);
    ket_same_ = (original_bs3_ == original_bs4_);
    braket_same_ = (original_bs1_ == original_bs3_ && original_bs2_ == original_bs4_);

    // Setup sieve data
    screening_threshold_ = Process::environment.options.get_double("INTS_TOLERANCE");
    auto screentype = Process::environment.options.get_str("SCREENING");
    if (screentype == "SCHWARZ")
        screening_type_ = ScreeningType::Schwarz;
    else if (screentype == "CSAM")
        screening_type_ = ScreeningType::CSAM;
    else if (screentype == "DENSITY")
        screening_type_ = ScreeningType::Density;
    else if (screentype == "NONE")
        screening_type_ = ScreeningType::None;
    else
        throw PSIEXCEPTION("Unknown screening type " + screentype + " in TwoBodyAOInt()");
    
    if (screening_threshold_ == 0.0) screening_type_ = ScreeningType::Schwarz;

}

TwoBodyAOInt::TwoBodyAOInt(const TwoBodyAOInt &rhs) : TwoBodyAOInt(rhs.integral_, rhs.deriv_) {
    buffers_.resize(rhs.buffers_.size());
    blocks12_ = rhs.blocks12_;
    blocks34_ = rhs.blocks34_;
    bra_same_ = rhs.bra_same_;
    ket_same_ = rhs.ket_same_;
    braket_same_ = rhs.braket_same_;
    screening_threshold_ = rhs.screening_threshold_;
    screening_threshold_squared_ = rhs.screening_threshold_squared_;
    nshell_ = rhs.nshell_;
    nbf_ = rhs.nbf_;
    screening_type_ = rhs.screening_type_;
    function_pair_values_ = rhs.function_pair_values_;
    shell_pair_values_ = rhs.shell_pair_values_;
    max_dens_shell_pair_ = rhs.max_dens_shell_pair_;
    shell_pair_exchange_values_ = rhs.shell_pair_exchange_values_;
    function_sqrt_ = rhs.function_sqrt_;
    function_pairs_ = rhs.function_pairs_;
    shell_pairs_ = rhs.shell_pairs_;
    shell_pairs_bra_ = rhs.shell_pairs_bra_;
    shell_pairs_ket_ = rhs.shell_pairs_ket_;
    max_integral_ = rhs.max_integral_;
    function_pairs_reverse_ = rhs.function_pairs_reverse_;
    shell_pairs_reverse_ = rhs.shell_pairs_reverse_;
    shell_to_shell_ = rhs.shell_to_shell_;
    function_to_function_ = rhs.function_to_function_;
    sieve_impl_ = rhs.sieve_impl_;
}

TwoBodyAOInt::~TwoBodyAOInt() {}

// Haser 1989, Equation 7 
void TwoBodyAOInt::update_density(const std::vector<SharedMatrix>& D) {

    if (max_dens_shell_pair_.size() == 0) {
        max_dens_shell_pair_.resize(D.size());
        for (int i = 0; i < D.size(); i++) {
            max_dens_shell_pair_[i].resize(nshell_ * nshell_);
        }
    }
    
    timer_on("Update Density");
#pragma omp parallel for
    for (int M = 0; M < nshell_; M++) {
        for (int N = M; N < nshell_; N++) {
            int m_start = bs1_->shell(M).function_index();
            int num_m = bs1_->shell(M).nfunction();

            int n_start = bs1_->shell(N).function_index();
            int num_n = bs1_->shell(N).nfunction();

            for (int i = 0; i < D.size(); i++) {
                double** Dp = D[i]->pointer();
                double max_dens = 0.0;
                for (int m = m_start; m < m_start + num_m; m++) {
                    for (int n = n_start; n < n_start + num_n; n++) {
                        max_dens = std::max(max_dens, std::abs(Dp[m][n]));
                    }
                }
                max_dens_shell_pair_[i][M * nshell_ + N] = max_dens;
                if (M != N) max_dens_shell_pair_[i][N * nshell_ + M] = max_dens;
            }

        }
    }
    timer_off("Update Density");

}


double TwoBodyAOInt::shell_pair_max_density(int M, int N) const {
    if (max_dens_shell_pair_.empty()) {
        throw PSIEXCEPTION("The density matrix has not been set in the TwoBodyAOInt class!");
    }
    double D_max = 0.0;
    for (const auto& matrix_max_per_pair: max_dens_shell_pair_) {
        D_max = std::max(D_max, matrix_max_per_pair[M * nshell_ + N]);
    }
    return D_max;
}

// Haser 1989 Equations 6 to 14
bool TwoBodyAOInt::shell_significant_density(int M, int N, int R, int S) {

    // Maximum density matrix equation
    double max_density = 0.0;

    // Equation 6 (RHF Case)
    if (max_dens_shell_pair_.size() == 1) {
        max_density = std::max({4.0 * max_dens_shell_pair_[0][M * nshell_ + N], 4.0 * max_dens_shell_pair_[0][R * nshell_ + S], 
            max_dens_shell_pair_[0][M * nshell_ + R], max_dens_shell_pair_[0][M * nshell_ + S],
            max_dens_shell_pair_[0][N * nshell_ + R], max_dens_shell_pair_[0][N * nshell_ + S]});
    } else { // UHF and ROHF
        // J-like terms
        double D_MN = max_dens_shell_pair_[0][M * nshell_ + N] + max_dens_shell_pair_[1][M * nshell_ + N];
        double D_RS = max_dens_shell_pair_[0][R * nshell_ + S] + max_dens_shell_pair_[1][R * nshell_ + S];

        // K-like terms
        double D_MR = std::max(max_dens_shell_pair_[0][M * nshell_ + R], max_dens_shell_pair_[1][M * nshell_ + R]);
        double D_MS = std::max(max_dens_shell_pair_[0][M * nshell_ + S], max_dens_shell_pair_[1][M * nshell_ + S]);
        double D_NR = std::max(max_dens_shell_pair_[0][N * nshell_ + R], max_dens_shell_pair_[1][N * nshell_ + R]);
        double D_NS = std::max(max_dens_shell_pair_[0][N * nshell_ + S], max_dens_shell_pair_[1][N * nshell_ + S]);

        max_density = std::max({2.0 * D_MN, 2.0 * D_RS, D_MR, D_MS, D_NR, D_NS});
    }

    // Square of Cauchy-Schwarz Q_MN terms (Eq. 13)
    double mn_mn = shell_pair_values_[N * nshell_ + M];
    double rs_rs = shell_pair_values_[S * nshell_ + R];

    // The density screened ERI bound (Eq. 6)
    return (mn_mn * rs_rs * max_density * max_density >= screening_threshold_squared_);
}

bool TwoBodyAOInt::shell_significant_csam(int M, int N, int R, int S) { 
    // Square of standard Cauchy-Schwarz Q_mu_nu terms (Eq. 1)
    double mn_mn = shell_pair_values_[N * nshell_ + M];
    double rs_rs = shell_pair_values_[S * nshell_ + R];

    // Square of M~_mu_nu terms (Eq. 9)
    double mm_rr = shell_pair_exchange_values_[R * nshell_ + M];
    double nn_ss = shell_pair_exchange_values_[S * nshell_ + N];
    double mm_ss = shell_pair_exchange_values_[S * nshell_ + M];
    double nn_rr = shell_pair_exchange_values_[R * nshell_ + N];

    // Square of M_mu_nu_lam_sig (Eq. 12)
    double csam_2 = std::max(mm_rr * nn_ss, mm_ss * nn_rr);

    // Square of Eq. 11
    double mnrs_2 = mn_mn * rs_rs * csam_2;

    return std::abs(mnrs_2) >= screening_threshold_squared_;
}

bool TwoBodyAOInt::shell_significant_schwarz(int M, int N, int R, int S) {
    return shell_pair_values_[N * nshell_ + M] * shell_pair_values_[R * nshell_ + S] >= screening_threshold_squared_;
}
bool TwoBodyAOInt::shell_significant_none(int M, int N, int R, int S) { return true; }

void TwoBodyAOInt::setup_sieve() {
    /*
     * Create metadata to screen (PQ|RS) using some relationship involving (PQ|PQ) and (RS|RS). We make a
     * few assumptions here: 1) we only want to create sieves where there is an equivalent basis set pair
     * in bra or ket.  So DF integrals like (A0|RS) with 0 denoting the dummy index will have a sieve
     * constructed from the RS pair in the ket that can be used to build dense indexing of the RS target.
     * The other type of DF integral (A0|B0) will not have any sieve associated with it.  Where dummy indices
     * are present they must occupy the 2nd and/or 4th index.
     */

    // Different approaches are possible, so we use a function pointer and set it once, to avoid logic later on
    switch (screening_type_) {
        case ScreeningType::CSAM:
            sieve_impl_ = [this](int M, int N, int R, int S) { return this->shell_significant_csam(M, N, R, S); };
            break;
        case ScreeningType::Schwarz:
            sieve_impl_ = [this](int M, int N, int R, int S) { return this->shell_significant_schwarz(M, N, R, S); };
            break;
        case ScreeningType::Density:
            sieve_impl_ = [this](int M, int N, int R, int S) { return this->shell_significant_density(M, N, R, S); };
            break;
        case ScreeningType::None:   
            sieve_impl_ = [this](int M, int N, int R, int S) { return this->shell_significant_none(M, N, R, S); };
            return;
        default:
            throw PSIEXCEPTION("Unimplemented screening type in TwoBodyAOInt::setup_sieve()");
    }


    // We assume that only the bra or the ket has a pair that generates a sieve.  If all bases are the same, either
    // can be used.  If only bra or ket has a matching pair, that matching pair is used.  If both bra and ket have
    // matching pairs but those pairs are different, we need to generalize this machinery a little to disambiguate
    // which pair should be used to form the sieve.  I don't know of a need for that right now, so I'll assume its
    // not needed and add a safety check to futureproof the code against that kind of use case further down the road.
    if(bra_same_ && ket_same_ && !braket_same_) throw PSIEXCEPTION("Unexpected integral type (aa|bb) in setup_sieve()");

    if(bra_same_) {
        create_sieve_pair_info(basis1(), shell_pairs_bra_, true);
        shell_pairs_ = shell_pairs_bra_;
    } else {
        if (basis2()->l2_shell(0) != libint2::Shell::unit()) 
               throw PSIEXCEPTION("If different basis sets exist in the bra, basis3 is expected to be dummy in setup_sieve()");
        for(int shell = 0; shell < basis1()->nshell(); ++shell) shell_pairs_bra_.emplace_back(shell,0);
    }
    if(ket_same_) {
        if(braket_same_) {
            shell_pairs_ket_ = shell_pairs_bra_;
        } else {
            create_sieve_pair_info(basis3(), shell_pairs_ket_, false);
            shell_pairs_ = shell_pairs_ket_;
        }
    } else {
        if (basis4()->l2_shell(0) != libint2::Shell::unit()) 
               throw PSIEXCEPTION("If different basis sets exist in the ket, basis4 is expected to be dummy in setup_sieve()");
        for(int shell = 0; shell < basis3()->nshell(); ++shell) shell_pairs_ket_.emplace_back(shell,0);
    }
}

void TwoBodyAOInt::create_sieve_pair_info(const std::shared_ptr<BasisSet> bs, PairList &shell_pairs, bool is_bra) {

    nshell_ = bs->nshell();
    nbf_ = bs->nbf();

    function_pair_values_.resize((size_t)nbf_ * nbf_, 0.0);
    shell_pair_values_.resize((size_t)nshell_ * nshell_, 0.0);
    max_integral_ = 0.0;

    bs1_ = bs;
    bs2_ = bs;
    bs3_ = bs;
    bs4_ = bs;
    for (int P = 0; P < nshell_; P++) {
        for (int Q = 0; Q <= P; Q++) {
            int nP = bs->shell(P).nfunction();
            int nQ = bs->shell(Q).nfunction();
            int oP = bs->shell(P).function_index();
            int oQ = bs->shell(Q).function_index();
            compute_shell_for_sieve(bs, P, Q, P, Q, is_bra);
            const double *buffer = target_full_;
            double shell_max_val = 0.0;
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    shell_max_val =
                        std::max(shell_max_val, std::abs(buffer[p * (nQ * nP * nQ + nQ) + q * (nP * nQ + 1)]));
                }
            }
            max_integral_ = std::max(max_integral_, shell_max_val);
            shell_pair_values_[P * nshell_ + Q] = shell_pair_values_[Q * nshell_ + P] = shell_max_val;
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    function_pair_values_[(p + oP) * nbf_ + (q + oQ)] = function_pair_values_[(q + oQ) * nbf_ + (p + oP)] = shell_max_val;
                }
            }
        }
    }
    bs1_ = original_bs1_;
    bs2_ = original_bs2_;
    bs3_ = original_bs3_;
    bs4_ = original_bs4_;

    screening_threshold_squared_ = screening_threshold_ * screening_threshold_;
    double screening_threshold_over_max = screening_threshold_ / max_integral_;
    double screening_threshold_squared_over_max = screening_threshold_squared_ / max_integral_;

    shell_pairs.clear();
    function_pairs_.clear();
    shell_pairs_reverse_.resize(nshell_ * (nshell_ + 1L) / 2L);
    function_pairs_reverse_.resize(nbf_ * (nbf_ + 1L) / 2L);

    long int offset = 0L;
    size_t munu = 0L;
    for (int mu = 0; mu < nbf_; mu++) {
        for (int nu = 0; nu <= mu; nu++, munu++) {
            if (function_pair_values_[mu * nbf_ + nu] >= screening_threshold_squared_over_max) {
                function_pairs_.push_back(std::make_pair(mu, nu));
                function_pairs_reverse_[munu] = offset;
                offset++;
            } else {
                function_pairs_reverse_[munu] = -1L;
            }
        }
    }

    shell_to_shell_.clear();
    function_to_function_.clear();
    shell_to_shell_.resize(nshell_);
    function_to_function_.resize(nbf_);

    for (int MU = 0; MU < nshell_; MU++) {
        for (int NU = 0; NU < nshell_; NU++) {
            if (shell_pair_values_[MU * nshell_ + NU] >= screening_threshold_squared_over_max) {
                shell_to_shell_[MU].push_back(NU);
            }
        }
    }

    shell_pairs.clear();
    std::fill_n(shell_pairs_reverse_.begin(), nshell_ * (nshell_ + 1) / 2, -1);

    offset = 0L;
    size_t MUNU = 0L;
    for (int MU = 0; MU < nshell_; MU++) {
        for (int NU = 0; NU <= MU; NU++, MUNU++) {
            if (shell_pair_values_[MU * nshell_ + NU] >= screening_threshold_squared_over_max) {
                shell_pairs.push_back(std::make_pair(MU, NU));
                shell_pairs_reverse_[MUNU] = offset;
                offset++;
            }
        }
    }

    for (int mu = 0; mu < nbf_; mu++) {
        for (int nu = 0; nu < nbf_; nu++) {
            if (function_pair_values_[mu * nbf_ + nu] >= screening_threshold_squared_over_max) {
                function_to_function_[mu].push_back(nu);
            }
        }
    }

    if (screening_type_ == ScreeningType::CSAM) {
        // Setup information for exchange term screening
        function_sqrt_.resize(nbf_);
        shell_pair_exchange_values_.resize((size_t)nshell_ * nshell_);
        std::fill(function_sqrt_.begin(), function_sqrt_.end(), 0.0);
        std::fill(shell_pair_exchange_values_.begin(), shell_pair_exchange_values_.end(), 0.0);

        for (int P = 0; P < nshell_; P++) {
            for (int Q = P; Q >= 0; Q--) {
                int nP = bs->shell(P).nfunction();
                int nQ = bs->shell(Q).nfunction();
                int oP = bs->shell(P).function_index();
                int oQ = bs->shell(Q).function_index();
                compute_shell_for_sieve(bs, P, P, Q, Q, is_bra);
                const double *buffer = target_full_;

                // Computing Q_mu_mu (for denominator of Eq.9)
                if (Q == P) {
                    int oP = bs->shell(P).function_index();
                    for (int p = 0; p < nP; ++p) {
                        function_sqrt_[oP + p] = std::sqrt(std::abs(buffer[p * (nP * nP * nP + nP) + p * (nP * nP + 1)]));
                    }
                }

                // Computing square of M~_mu_lam (defined in Eq. 9)
                double max_val = 0.0;
                for (int p = 0; p < nP; p++) {
                    for (int q = 0; q < nQ; q++) {
                        max_val = std::max(max_val, std::abs(buffer[p * nQ * nQ * (nP + 1) + q * (nQ + 1)]) /
                                                        (function_sqrt_[p + oP] * function_sqrt_[q + oQ]));
                    }
                }
                shell_pair_exchange_values_[P * nshell_ + Q] = shell_pair_exchange_values_[Q * nshell_ + P] = max_val;
            }
        }
    }
}

std::shared_ptr<BasisSet> TwoBodyAOInt::basis() { return original_bs1_; }

std::shared_ptr<BasisSet> TwoBodyAOInt::basis1() { return original_bs1_; }

std::shared_ptr<BasisSet> TwoBodyAOInt::basis2() { return original_bs2_; }

std::shared_ptr<BasisSet> TwoBodyAOInt::basis3() { return original_bs3_; }

std::shared_ptr<BasisSet> TwoBodyAOInt::basis4() { return original_bs4_; }

std::vector<ShellPairBlock> TwoBodyAOInt::get_blocks12() const { return blocks12_; }

std::vector<ShellPairBlock> TwoBodyAOInt::get_blocks34() const { return blocks34_; }

// Expected to be overridden by derived classes
void TwoBodyAOInt::create_blocks() {
    // Default implementation : do no blocking but use the sieved shell pairs if possible
    // Each ShellPairBlock will only contain one shell pair
    blocks12_.clear();
    blocks34_.clear();

    const auto nshell1 = basis1()->nshell();
    const auto nshell2 = basis2()->nshell();
    const auto nshell3 = basis3()->nshell();
    const auto nshell4 = basis4()->nshell();

    // Push back only the pairs that survived the sieving process.  This is only
    // possible if all four basis sets are the same in the current implementation.
    blocks12_.reserve(shell_pairs_bra_.size());
    for (const auto &pair : shell_pairs_bra_) {
        const auto &s1 = pair.first;
        const auto &s2 = pair.second;
        blocks12_.push_back({{s1, s2}});
    }
    blocks34_.reserve(shell_pairs_ket_.size());
    for (const auto &pair : shell_pairs_ket_) {
        const auto &s3 = pair.first;
        const auto &s4 = pair.second;
        blocks34_.push_back({{s3, s4}});
    }
}

bool TwoBodyAOInt::shell_block_significant(int shellpair12, int shellpair34) const {
    const auto &vsh12 = blocks12_[shellpair12];
    const auto &vsh34 = blocks34_[shellpair34];

    for (const auto &sh12 : vsh12) {
        for (const auto &sh34 : vsh34) {
            if(shell_significant(sh12.first, sh12.second, sh34.first, sh34.second)) return true;
        }
    }

    return false;
}

bool TwoBodyAOInt::shell_pair_significant(int M, int N) const {
    return shell_pair_values_[M * nshell_ + N] * max_integral_ >= screening_threshold_squared_;
}

void TwoBodyAOInt::compute_shell_blocks(int shellpair12, int shellpair34, int npair12, int npair34) {
    // Default implementation - go through the blocks and do each quartet
    // one at a time

    // reset the target & source pointers
    target_ = target_full_;
    source_ = source_full_;

    const auto &vsh12 = blocks12_[shellpair12];
    const auto &vsh34 = blocks34_[shellpair34];

    for (const auto &sh12 : vsh12) {
        const auto &shell1 = bs1_->shell(sh12.first);
        const auto &shell2 = bs2_->shell(sh12.second);

        int n1 = shell1.nfunction();
        int n2 = shell2.nfunction();

        for (const auto &sh34 : vsh34) {
            const auto &shell3 = bs3_->shell(sh34.first);
            const auto &shell4 = bs4_->shell(sh34.second);

            int n3 = shell3.nfunction();
            int n4 = shell4.nfunction();

            const int n1234 = n1 * n2 * n3 * n4;

            // actually compute the eri
            // this will put the results in target_
            auto ret = compute_shell(sh12.first, sh12.second, sh34.first, sh34.second);

            //// advance the target pointer
            target_ += n1234;

            // Since we are only doing one at a time we don't need to
            // move the source_ pointer
        }
    }
}

void TwoBodyAOInt::compute_shell_blocks_deriv1(int shellpair12, int shellpair34, int npair12, int npair34) {
    // Default implementation - go through the blocks and do each quartet
    // one at a time

    // reset the target & source pointers
    target_ = target_full_;
    source_ = source_full_;

    auto vsh12 = blocks12_[shellpair12];
    auto vsh34 = blocks34_[shellpair34];

    for (const auto sh12 : vsh12) {
        const auto &shell1 = original_bs1_->shell(sh12.first);
        const auto &shell2 = original_bs2_->shell(sh12.second);

        int n1 = shell1.nfunction();
        int n2 = shell2.nfunction();

        for (const auto sh34 : vsh34) {
            const auto &shell3 = original_bs3_->shell(sh34.first);
            const auto &shell4 = original_bs4_->shell(sh34.second);

            int n3 = shell3.nfunction();
            int n4 = shell4.nfunction();

            const int n1234 = 12 * n1 * n2 * n3 * n4;

            // actually compute the eri
            // this will put the results in target_
            compute_shell_deriv1(sh12.first, sh12.second, sh34.first, sh34.second);

            //// advance the target pointer
            target_ += n1234;

            // Since we are only doing one at a time we don't need to
            // move the source_ pointer
        }
    }
}

void TwoBodyAOInt::compute_shell_blocks_deriv2(int shellpair12, int shellpair34, int npair12, int npair34) {
    // Default implementation - go through the blocks and do each quartet
    // one at a time

    // reset the target & source pointers
    target_ = target_full_;
    source_ = source_full_;

    auto vsh12 = blocks12_[shellpair12];
    auto vsh34 = blocks34_[shellpair34];

    for (const auto sh12 : vsh12) {
        const auto &shell1 = original_bs1_->shell(sh12.first);
        const auto &shell2 = original_bs2_->shell(sh12.second);

        int n1 = shell1.nfunction();
        int n2 = shell2.nfunction();

        for (const auto sh34 : vsh34) {
            const auto &shell3 = original_bs3_->shell(sh34.first);
            const auto &shell4 = original_bs4_->shell(sh34.second);

            int n3 = shell3.nfunction();
            int n4 = shell4.nfunction();

            // TODO: figure out if tranlational invariance relations are applied automagically
            const int n1234 = 78 * n1 * n2 * n3 * n4;

            // actually compute the eri
            // this will put the results in target_
            compute_shell_deriv2(sh12.first, sh12.second, sh34.first, sh34.second);

            //// advance the target pointer
            target_ += n1234;

            // Since we are only doing one at a time we don't need to
            // move the source_ pointer
        }
    }
}

void TwoBodyAOInt::normalize_am(std::shared_ptr<GaussianShell> s1, std::shared_ptr<GaussianShell> s2,
                                std::shared_ptr<GaussianShell> s3, std::shared_ptr<GaussianShell> s4, int nchunk) {
    // Integrals assume this normalization is 1.0.
    return;
#if 0
#ifdef MINTS_TIMER
    timer_on("Angular momentum normalization");
#endif
    int am1 = s1->am();
    int am2 = s2->am();
    int am3 = s3->am();
    int am4 = s4->am();
    int length = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

    // Need to go through and grab all the integrals for this given shell and add them
    // to the running totals.
    int nprim = 0;
    for (int i = 0; i <= am1; i++) {
        int l1 = am1 - i;
        for (int j = 0; j <= i; j++) {
            int m1 = i - j;
            int n1 = j;
            double norm_a = s1->normalize(l1, m1, n1);

            for (int k = 0; k <= am2; k++) {
                int l2 = am2 - k;
                for (int l = 0; l <= k; l++) {
                    int m2 = k - l;
                    int n2 = l;
                    double norm_b = s2->normalize(l2, m2, n2);

                    for (int m = 0; m <= am3; m++) {
                        int l3 = am3 - m;
                        for (int n = 0; n <= m; n++) {
                            int m3 = m - n;
                            int n3 = n;
                            double norm_c = s3->normalize(l3, m3, n3);

                            for (int o = 0; o <= am4; o++) {
                                int l4 = am4 - o;
                                for (int p = 0; p <= o; p++) {
                                    int m4 = o - p;
                                    int n4 = p;
                                    double norm_d = s4->normalize(l4, m4, n4);

                                    // printf("normalization %f %f %f %f\n", norm_a, norm_b, norm_c, norm_d);
                                    for (int chunk=0; chunk < nchunk; ++chunk) {
                                        source_[nprim+(chunk*length)] *= norm_a * norm_b * norm_c * norm_d;
                                    }
                                    nprim++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("Angular momentum normalization");
#endif
#endif
}

void TwoBodyAOInt::permute_target(double *s, double *t, int sh1, int sh2, int sh3, int sh4, bool p12, bool p34,
                                  bool p13p24) {
#ifdef MINTS_TIMER
    timer_on("Permute target");
#endif
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);
    const GaussianShell &s3 = bs3_->shell(sh3);
    const GaussianShell &s4 = bs4_->shell(sh4);

    int nbf1 = s1.nfunction();
    int nbf2 = s2.nfunction();
    int nbf3 = s3.nfunction();
    int nbf4 = s4.nfunction();

    if (!p13p24) {
        if (p12) {
            if (p34) {
                permute_1234_to_2143(s, t, nbf1, nbf2, nbf3, nbf4);
            } else {
                permute_1234_to_2134(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        } else {
            permute_1234_to_1243(s, t, nbf1, nbf2, nbf3, nbf4);
        }
    } else {
        if (p12) {
            if (p34) {
                permute_1234_to_4321(s, t, nbf1, nbf2, nbf3, nbf4);
            } else {
                permute_1234_to_4312(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        } else {
            if (p34) {
                permute_1234_to_3421(s, t, nbf1, nbf2, nbf3, nbf4);
            } else {
                permute_1234_to_3412(s, t, nbf1, nbf2, nbf3, nbf4);
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("Permute target");
#endif
}

void TwoBodyAOInt::permute_1234_to_1243(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf1;
    int f2 = nbf2;
    int f3 = nbf4;
    int f4 = nbf3;
    for (int bf1 = 0; bf1 < f1; bf1++) {
        for (int bf2 = 0; bf2 < f2; bf2++) {
            for (int bf4 = 0; bf4 < f4; bf4++) {
                for (int bf3 = 0; bf3 < f3; bf3++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_2134(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf2;
    int f2 = nbf1;
    int f3 = nbf3;
    int f4 = nbf4;
    for (int bf2 = 0; bf2 < f2; bf2++) {
        for (int bf1 = 0; bf1 < f1; bf1++) {
            for (int bf3 = 0; bf3 < f3; bf3++) {
                for (int bf4 = 0; bf4 < f4; bf4++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_2143(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf2;
    int f2 = nbf1;
    int f3 = nbf4;
    int f4 = nbf3;
    for (int bf2 = 0; bf2 < f2; bf2++) {
        for (int bf1 = 0; bf1 < f1; bf1++) {
            for (int bf4 = 0; bf4 < f4; bf4++) {
                for (int bf3 = 0; bf3 < f3; bf3++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_3412(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf3;
    int f2 = nbf4;
    int f3 = nbf1;
    int f4 = nbf2;
    for (int bf3 = 0; bf3 < f3; bf3++) {
        for (int bf4 = 0; bf4 < f4; bf4++) {
            for (int bf1 = 0; bf1 < f1; bf1++) {
                for (int bf2 = 0; bf2 < f2; bf2++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_4312(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf4;
    int f2 = nbf3;
    int f3 = nbf1;
    int f4 = nbf2;
    for (int bf3 = 0; bf3 < f3; bf3++) {
        for (int bf4 = 0; bf4 < f4; bf4++) {
            for (int bf2 = 0; bf2 < f2; bf2++) {
                for (int bf1 = 0; bf1 < f1; bf1++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_3421(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf3;
    int f2 = nbf4;
    int f3 = nbf2;
    int f4 = nbf1;
    for (int bf4 = 0; bf4 < f4; bf4++) {
        for (int bf3 = 0; bf3 < f3; bf3++) {
            for (int bf1 = 0; bf1 < f1; bf1++) {
                for (int bf2 = 0; bf2 < f2; bf2++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::permute_1234_to_4321(double *s, double *t, int nbf1, int nbf2, int nbf3, int nbf4) {
    int f1 = nbf4;
    int f2 = nbf3;
    int f3 = nbf2;
    int f4 = nbf1;
    for (int bf4 = 0; bf4 < f4; bf4++) {
        for (int bf3 = 0; bf3 < f3; bf3++) {
            for (int bf2 = 0; bf2 < f2; bf2++) {
                for (int bf1 = 0; bf1 < f1; bf1++) {
                    double *t_ptr = t + ((bf1 * f2 + bf2) * f3 + bf3) * f4 + bf4;
                    *(t_ptr) = *(s++);
                }
            }
        }
    }
}

void TwoBodyAOInt::pure_transform(int sh1, int sh2, int sh3, int sh4, int nchunk, bool copy_to_source) {
#ifdef MINTS_TIMER
    timer_on("Pure transformation");
#endif
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);
    const GaussianShell &s3 = bs3_->shell(sh3);
    const GaussianShell &s4 = bs4_->shell(sh4);

    // Get the transforms from the basis set
    SphericalTransformIter trans1(*integral()->spherical_transform(s1.am()));
    SphericalTransformIter trans2(*integral()->spherical_transform(s2.am()));
    SphericalTransformIter trans3(*integral()->spherical_transform(s3.am()));
    SphericalTransformIter trans4(*integral()->spherical_transform(s4.am()));

    // Get the angular momentum for each shell
    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();

    // Get number of Cartesian functions for each shell
    int nao1 = s1.ncartesian();
    int nao2 = s2.ncartesian();
    int nao3 = s3.ncartesian();
    int nao4 = s4.ncartesian();

    int nbf1 = s1.nfunction();
    int nbf2 = s2.nfunction();
    int nbf3 = s3.nfunction();
    int nbf4 = s4.nfunction();

    // Get if each shell has pure functions
    bool is_pure1 = s1.is_pure();
    bool is_pure2 = s2.is_pure();
    bool is_pure3 = s3.is_pure();
    bool is_pure4 = s4.is_pure();

    for (int ichunk = 0; ichunk < nchunk; ++ichunk) {
        // Compute the offset in source_, and target
        size_t sourcechunkoffset = ichunk * (nao1 * nao2 * nao3 * nao4);
        size_t targetchunkoffset = ichunk * ((size_t)nbf1 * nbf2 * nbf3 * nbf4);
        double *source1, *target1;
        double *source2, *target2;
        double *source3, *target3;
        double *source4, *target4;
        double *source = source_ + sourcechunkoffset;
        double *target = target_ + targetchunkoffset;
        double *tmpbuf = tformbuf_;

        int transform_index = 8 * is_pure1 + 4 * is_pure2 + 2 * is_pure3 + is_pure4;
        switch (transform_index) {
            case 0:
                break;

            case 1:
                source4 = source;
                target4 = target;
                break;

            case 2:
                source3 = source;
                target3 = target;
                break;

            case 3:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = target;
                break;

            case 4:
                source2 = source;
                target2 = target;
                break;

            case 5:
                source4 = source;
                target4 = tmpbuf;
                source2 = tmpbuf;
                target2 = target;
                break;

            case 6:
                source3 = source;
                target3 = tmpbuf;
                source2 = tmpbuf;
                target2 = target;
                break;

            case 7:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = source;
                source2 = source;
                target2 = target;
                break;

            case 8:
                source1 = source;
                target1 = target;
                break;

            case 9:
                source4 = source;
                target4 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;

            case 10:
                source3 = source;
                target3 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;

            case 11:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = source;
                source1 = source;
                target1 = target;
                break;

            case 12:
                source2 = source;
                target2 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;

            case 13:
                source4 = source;
                target4 = tmpbuf;
                source2 = tmpbuf;
                target2 = source;
                source1 = source;
                target1 = target;
                break;

            case 14:
                source3 = source;
                target3 = tmpbuf;
                source2 = tmpbuf;
                target2 = source;
                source1 = source;
                target1 = target;
                break;

            case 15:
                source4 = source;
                target4 = tmpbuf;
                source3 = tmpbuf;
                target3 = source;
                source2 = source;
                target2 = tmpbuf;
                source1 = tmpbuf;
                target1 = target;
                break;
        }

        size_t size = nbf1 * nbf2 * nbf3 * nbf4;
        if (is_pure4) {
            transform2e_4(am4, trans4, source4, target4, nao1 * nao2 * nao3, nao4);
            //            size *= nbf4;
        }
        if (is_pure3) {
            transform2e_3(am3, trans3, source3, target3, nao1 * nao2, nao3, nbf4);
            //            size *= nbf3;
        }
        if (is_pure2) {
            transform2e_2(am2, trans2, source2, target2, nao1, nao2, nbf3 * nbf4);
            //            size *= nbf2;
        }
        if (is_pure1) {
            transform2e_1(am1, trans1, source1, target1, nbf2 * nbf3 * nbf4);
            //            size *= nbf1;
        }
        // The permute indices routines depend on the integrals being in source_
        if (copy_to_source && (is_pure1 || is_pure2 || is_pure3 || is_pure4))
            memcpy(source, target, size * sizeof(double));
    }
#ifdef MINTS_TIMER
    timer_off("Pure transformation");
#endif
}

static void transform2e_1(int am, SphericalTransformIter &sti, double *s, double *t, int njkl) {
    memset(t, 0, INT_NPURE(am) * njkl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex() * njkl;
        double *tptr = t + sti.pureindex() * njkl;
        double coef = sti.coef();

        //        outfile->Printf( "2e_1: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(),
        //        sti.coef());

        for (int jkl = 0; jkl < njkl; jkl++) *(tptr++) += coef * *(sptr++);
    }
}

static void transform2e_2(int am, SphericalTransformIter &sti, double *s, double *t, int ni, int nj, int nkl) {
    int sj = INT_NPURE(am);
    const int sjkl = nj * nkl;
    const int tjkl = sj * nkl;

    memset(t, 0, ni * tjkl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex() * nkl;
        double *tptr = t + sti.pureindex() * nkl;
        double coef = sti.coef();

        //        outfile->Printf( "2e_2: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(),
        //        sti.coef());

        for (int i = 0; i < ni; i++, sptr += sjkl, tptr += tjkl) {
            for (int kl = 0; kl < nkl; kl++) tptr[kl] += coef * sptr[kl];
        }
    }
}

static void transform2e_3(int am, SphericalTransformIter &sti, double *s, double *t, int nij, int nk, int nl) {
    int sk = INT_NPURE(am);
    const int skl = nk * nl;
    const int tkl = sk * nl;

    memset(t, 0, nij * tkl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        double *sptr = s + sti.cartindex() * nl;
        double *tptr = t + sti.pureindex() * nl;
        // printf("cartindex = %d, pureindex = %d\n", sti.cartindex(), sti.pureindex());

        //        outfile->Printf( "2e_3: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(),
        //        sti.coef());

        double coef = sti.coef();
        for (int ij = 0; ij < nij; ij++, sptr += skl, tptr += tkl) {
            for (int l = 0; l < nl; l++) tptr[l] += coef * sptr[l];
        }
    }
}

// am => angular momentum of l
// sti => spherical transformation iterator
// s => source integrals buffer
// t => target buffer
// nijk => how many i, j, k combinations are there?
// nl => how man l's are there?
static void transform2e_4(int am, SphericalTransformIter &sti, double *s, double *t, int nijk, int nl) {
    // Protect ourselves
    const int sl = nl;
    const int tl = INT_NPURE(am);

    // Clear out target memory
    memset(t, 0, nijk * tl * sizeof(double));

    for (sti.first(); !sti.is_done(); sti.next()) {
        // Starting point in source and target buffers
        double *sptr = s + sti.cartindex();
        double *tptr = t + sti.pureindex();

        //        outfile->Printf( "2e_4: cart = %d pure = %d coef = %8.5f\n", sti.cartindex(), sti.pureindex(),
        //        sti.coef());

        // What's the coefficient we're using
        double coef = sti.coef();
        for (int ijk = 0; ijk < nijk; ++ijk) {
            // Add contribution of the source to the target
            *(tptr) += coef * *(sptr);

            // skip ahead to the next ijk
            sptr += sl;
            tptr += tl;
        }
    }
}
