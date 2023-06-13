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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/fjt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <libint2/shell.h>
#include <libint2/engine.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

// libderiv computes 9 of the 12 total derivatives. It computes 3 of the
// centers we handle the 4th.
#define ERI_1DER_NTYPE (9)
// libderiv computes the second derivatives using the first derivatives.
// The first derivatives are provided when second derivatives are asked
// for.
#define ERI_2DER_NTYPE (ERI_1DER_NTYPE + 45)

using namespace psi;

namespace {

unsigned char ntypes[] = {1, ERI_1DER_NTYPE, ERI_2DER_NTYPE};

}  // end namespace

///// Libint2 implementation

Libint2TwoElectronInt::Libint2TwoElectronInt(const IntegralFactory *integral, int deriv, double screening_threshold,
                                             bool use_shell_pairs, bool needs_exchange)
    : TwoBodyAOInt(integral, deriv), use_shell_pairs_(use_shell_pairs) {
    // Initialize libint static data

    // Make sure there's enough space for the sieve generation.  This array is used to return an array of
    // zeros back to the caller if libint2 gave us nullptr, so the caller doesn't have to check.
    size_t sieve_size = std::max(basis1()->max_function_per_shell() * basis2()->max_function_per_shell() *
                                     basis1()->max_function_per_shell() * basis2()->max_function_per_shell(),
                                 basis3()->max_function_per_shell() * basis4()->max_function_per_shell() *
                                     basis3()->max_function_per_shell() * basis4()->max_function_per_shell());
    size_t size = std::max((size_t)basis1()->max_function_per_shell() * basis2()->max_function_per_shell() *
                               basis3()->max_function_per_shell() * basis4()->max_function_per_shell(),
                           sieve_size);
    zero_vec_ = std::vector<double>(size, 0.0);
}

Libint2TwoElectronInt::Libint2TwoElectronInt(const Libint2TwoElectronInt &rhs)
    : TwoBodyAOInt(rhs),
      schwarz_engine_(rhs.schwarz_engine_),
      braket_(rhs.braket_),
      use_shell_pairs_(rhs.use_shell_pairs_) {
    pairs12_ = rhs.pairs12_;
    pairs34_ = rhs.pairs34_;
    zero_vec_ = rhs.zero_vec_;
    for (const auto &e : rhs.engines_) engines_.emplace_back(e);
}

void Libint2TwoElectronInt::common_init() {
    bool dummy1 = basis1()->l2_shell(0) == libint2::Shell::unit();
    bool dummy2 = basis2()->l2_shell(0) == libint2::Shell::unit();
    bool dummy3 = basis3()->l2_shell(0) == libint2::Shell::unit();
    bool dummy4 = basis4()->l2_shell(0) == libint2::Shell::unit();

    if (!dummy1 && !dummy2 && !dummy3 && !dummy4) {
        braket_ = libint2::BraKet::xx_xx;
    } else if (!dummy1 && dummy2 && !dummy3 && !dummy4) {
        braket_ = libint2::BraKet::xs_xx;
    } else if (!dummy1 && !dummy2 && !dummy3 && dummy4) {
        braket_ = libint2::BraKet::xx_xs;
    } else if (!dummy1 && dummy2 && !dummy3 && dummy4) {
        braket_ = libint2::BraKet::xs_xs;
    } else {
        throw PSIEXCEPTION("Bad BraKet type in Libint2TwoElectronInt");
    }

    for (auto &engine : engines_) engine.set(braket_);

    int num_chunks;
    switch (deriv_) {
        case 0:
            num_chunks = 1;
            break;
        case 1:
            num_chunks = 12;
            break;
        case 2:
            num_chunks = 78;
            break;
        default:
            throw PSIEXCEPTION("Libint2 engine only supports up to second derivatives currently.");
    }
    buffers_.resize(num_chunks);

    target_full_ = const_cast<double *>(engines_[0].results()[0]);
    target_ = target_full_;

    // Make sure the engine can handle the type of integral used to build a sieve
    setup_sieve();
    // Reset the engine type back to the general case needed
    create_blocks();
    const auto max_engine_precision = std::numeric_limits<double>::epsilon() * screening_threshold_;

    size_t npairs = shell_pairs_bra_.size();
    pairs12_.resize(npairs);
    // #pragma omp parallel for
    for (int pair = 0; pair < npairs; ++pair) {
        auto s1 = shell_pairs_bra_[pair].first;
        auto s2 = shell_pairs_bra_[pair].second;
        pairs12_[pair] = std::make_shared<libint2::ShellPair>(basis1()->l2_shell(s1), basis2()->l2_shell(s2),
                                                              std::log(max_engine_precision));
    }
    npairs = shell_pairs_ket_.size();
    pairs34_.resize(npairs);
    // #pragma omp parallel for
    for (int pair = 0; pair < npairs; ++pair) {
        auto s3 = shell_pairs_ket_[pair].first;
        auto s4 = shell_pairs_ket_[pair].second;
        pairs34_[pair] = std::make_shared<libint2::ShellPair>(basis3()->l2_shell(s3), basis4()->l2_shell(s4),
                                                              std::log(max_engine_precision));
    }
}

Libint2TwoElectronInt::~Libint2TwoElectronInt() {}

void Libint2TwoElectronInt::initialize_sieve() { 
    create_sieve_pair_info();
    
    create_blocks();
    const auto max_engine_precision = std::numeric_limits<double>::epsilon() * screening_threshold_;

    // Reset the engine type back to the general case needed
    size_t npairs = shell_pairs_bra_.size();
    pairs12_.resize(npairs);
//#pragma omp parallel for
    for (int pair = 0; pair < npairs; ++pair) {
        auto s1 = shell_pairs_bra_[pair].first;
        auto s2 = shell_pairs_bra_[pair].second;
        pairs12_[pair] = std::make_shared<libint2::ShellPair>(basis1()->l2_shell(s1), basis2()->l2_shell(s2),
                                                              std::log(max_engine_precision));
    }
    npairs = shell_pairs_ket_.size();
    pairs34_.resize(npairs);
//#pragma omp parallel for
    for (int pair = 0; pair < npairs; ++pair) {
        auto s3 = shell_pairs_ket_[pair].first;
        auto s4 = shell_pairs_ket_[pair].second;
        pairs34_[pair] = std::make_shared<libint2::ShellPair>(basis3()->l2_shell(s3), basis4()->l2_shell(s4),
                                                              std::log(max_engine_precision));
    }
}

size_t Libint2TwoElectronInt::compute_shell(const AOShellCombinationsIterator &shellIter) {
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t Libint2TwoElectronInt::compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3,
                                                      int s4, bool is_bra) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_for_sieve");
#endif
    const auto &sh1 = bs->l2_shell(s1);
    const auto &sh2 = bs->l2_shell(s2);
    const auto &sh3 = bs->l2_shell(s3);
    const auto &sh4 = bs->l2_shell(s4);

    schwarz_engine_.compute(sh1, sh2, sh3, sh4);

    size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();
    buffers_[0] = target_full_ = const_cast<double *>(schwarz_engine_.results()[0]);
    if (target_full_ == nullptr) {
        // The caller will try to read the buffer if there isn't a check on the number of ints computed
        // so we point to a valid array of zeros here to prevent memory bugs in the calling routine.
        buffers_[0] = target_full_ = zero_vec_.data();
        ntot = 0;
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_for_sieve");
#endif
    return ntot;
}

size_t Libint2TwoElectronInt::compute_shell(int s1, int s2, int s3, int s4) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell");
#endif

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    libint2_wrapper0(sh1, sh2, sh3, sh4);

    size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();

    buffers_[0] = target_full_ = const_cast<double *>(engines_[0].results()[0]);
    if (target_full_ == nullptr) {
        // The caller will try to read the buffer if there isn't a check on the number of ints computed
        // so we point to a valid array of zeros here to prevent memory bugs in the calling routine.
        buffers_[0] = target_full_ = zero_vec_.data();
        ntot = 0;
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell");
#endif
    return ntot;
}

size_t Libint2TwoElectronInt::compute_shell_deriv1(int s1, int s2, int s3, int s4) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_deriv1");
#endif

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    libint2_wrapper1(sh1, sh2, sh3, sh4);

    size_t ntot = 0;
    bool none_computed = engines_[1].results()[0] == nullptr;
    if (none_computed) {
        for (int i = 0; i < 12; ++i) {
            buffers_[i] = zero_vec_.data();
        }
    } else {
        for (int i = 0; i < 12; ++i) {
            if (i < engines_[1].results().size()) {
                buffers_[i] = engines_[1].results()[i];
            } else {
                buffers_[i] = zero_vec_.data();
            }
        }
        ntot = 12 * sh1.size() * sh2.size() * sh3.size() * sh4.size();
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_deriv1");
#endif
    return ntot;
}

size_t Libint2TwoElectronInt::compute_shell_deriv2(int s1, int s2, int s3, int s4) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_deriv2");
#endif

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    libint2_wrapper2(sh1, sh2, sh3, sh4);

    size_t ntot = 0;
    bool none_computed = engines_[2].results()[0] == nullptr;
    if (none_computed) {
        for (int i = 0; i < 78; ++i) {
            buffers_[i] = zero_vec_.data();
        }
    } else {
        for (int i = 0; i < 78; ++i) {
            if (i < engines_[2].results().size()) {
                buffers_[i] = engines_[2].results()[i];
            } else {
                buffers_[i] = zero_vec_.data();
            }
        }
        ntot = 78 * sh1.size() * sh2.size() * sh3.size() * sh4.size();
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_deriv2");
#endif
    return ntot;
}

void Libint2TwoElectronInt::compute_shell_blocks(int shellpair12, int shellpair34, int npair12, int npair34) {
    if (npair12 != -1 || npair34 != -1)
        throw PSIEXCEPTION("npair12 and npair34 arguments are not supported by the Libint2 engine.");
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_blocks");
#endif
    // This engine doesn't block shells, so each "block" is just 1 shell
    int s1 = blocks12_[shellpair12][0].first;
    int s2 = blocks12_[shellpair12][0].second;
    int s3 = blocks34_[shellpair34][0].first;
    int s4 = blocks34_[shellpair34][0].second;

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    size_t ntot = 0;

    const auto *sp12 = pairs12_[shellpair12].get();
    const auto *sp34 = pairs34_[shellpair34].get();
    libint2_wrapper0(sh1, sh2, sh3, sh4, sp12, sp34);

    target_full_ = const_cast<double *>(engines_[0].results()[0]);
    if (target_full_) {
        buffers_[0] = engines_[0].results()[0];
        ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();
    } else {
        target_full_ = zero_vec_.data();
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_blocks");
#endif
}
