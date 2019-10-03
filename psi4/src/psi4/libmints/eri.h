/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

#include <unordered_map>
#include <libint2.hpp>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/twobody.h"
#include "psi4/psi4-dec.h"

namespace psi {

class BasisSet;
class TwoBodyAOInt;
class IntegralFactory;
class Fjt;
class AOShellCombinationsIterator;
class CorrelationFactor;

/**
 * \ingroup MINTS
 * Structure to hold precomputed gaussian product information
 */
struct PrimPair {
    //! x, y, z coordinate of gaussian
    double P[3];
    //! Distance between P and shell i center
    double PA[3];
    //! Distance between P and shell j center
    double PB[3];
    //! alphas for both centers
    double ai, aj;
    //! gamma (ai + aj)
    double gamma;
    //! Contraction coefficients
    double ci, cj;
    //! Overlap between primitives
    double overlap;
};

/**
 * \ingroup MINTS
 * Structure to hold precomputed shell pair information
 */
struct ShellPair {
    //! Shells for this information.
    int i, j;
    //! Distance between shell i and shell j centers
    double AB[3];
    //! Vector of significant primitive pairs that define this shell pair
    std::vector<PrimPair> nonzeroPrimPairs;
};

/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class TwoElectronInt : public TwoBodyAOInt {
   protected:
    //! Libint object.
    Libint1_t libint_;
    //! Libderiv object
    Libderiv_t libderiv_;

    //! Maximum cartesian class size.
    int max_cart_;

    //! Computes the fundamental
    Fjt* fjt_;

    //! Computes the ERIs between four shells.
    size_t compute_quartet(int, int, int, int);

    //! Computes the ERI derivatives between four shells.
    size_t compute_quartet_deriv1(int, int, int, int);

    //! Computes the ERI second derivative between four shells.
    size_t compute_quartet_deriv2(int, int, int, int);

    //! Form shell pair information. Must be smart enough to handle arbitrary basis sets
    void init_shell_pairs12();
    void init_shell_pairs34();

    //! Should we use shell pair information?
    bool use_shell_pairs_;

    //! Shell pair information
    std::shared_ptr<std::vector<std::vector<ShellPair>>> pairs12_, pairs34_;

    //! Evaluates how much memory (in doubles) is needed to store shell pair data
    size_t memory_to_store_shell_pairs(const std::shared_ptr<BasisSet>&, const std::shared_ptr<BasisSet>&);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;

   public:
    //! Constructor. Use an IntegralFactory to create this object.
    TwoElectronInt(const IntegralFactory* integral, int deriv = 0, bool use_shell_pairs = false);

    ~TwoElectronInt() override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&) override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(int, int, int, int) override;

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv1(int, int, int, int) override;

    /// Compute ERI second derivatives between 4 sheels. Result is stored in buffer.
    size_t compute_shell_deriv2(int, int, int, int) override;
};

/*! \ingroup MINTS
 *  \class Libint2TwoElectronInt
 *  \brief Capable of computing two-electron repulsion integrals.
 */
template <libint2::Operator op, libint2::BraKet bk>
class Libint2TwoElectronInt : public TwoBodyAOInt {
    using ShellPairData = std::vector<std::shared_ptr<libint2::ShellPair>>;

   protected:
    //! Libint2 engine
    libint2::Engine libint2_;
    std::vector<double> results_;

    //! Form shell pair information
    std::tuple<std::vector<ShellPairBlock>, ShellPairData> build_shell_pair_list(std::shared_ptr<BasisSet> bs1,
                                                                                 std::shared_ptr<BasisSet> bs2) {
        const auto nsh1 = bs1->nshell();
        const auto nsh2 = bs2->nshell();
        const auto bs1_equiv_bs2 = (bs1 == bs2);
        std::vector<ShellPairBlock> blocks;
        auto nthreads = Process::environment.get_n_threads();

        // construct the 2-electron repulsion integrals engine
        using libint2::Engine;
        std::vector<Engine> engines;
        engines.reserve(nthreads);
        engines.emplace_back(libint2::Operator::overlap, std::max(bs1->max_nprimitive(), bs2->max_nprimitive()),
                             std::max(bs1->max_am(), bs2->max_am()), 0);
        for (size_t i = 1; i != nthreads; ++i) {
            engines.push_back(engines[0]);
        }

        std::mutex mx;

#pragma omp parallel
        {
            int thread_id = 0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif
            auto& engine = engines[thread_id];
            const auto& buf = engine.results();

            // loop over permutationally-unique set of shells
            for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
                auto n1 = bs1->shell(s1).nfunction();

                auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
                for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
                    if (s12 % nthreads != thread_id) continue;

                    auto on_same_center = (bs1->shell(s1).center() == bs2->shell(s2).center());
                    bool significant = on_same_center;
                    if (!on_same_center) {
                        auto n2 = bs2->shell(s2).nfunction();
                        engines[thread_id].compute(bs1->l2_shell(s1), bs2->l2_shell(s2));
                        double normsq = std::inner_product(buf[0], buf[0] + n1 * n2, buf[0], 0.0);
                        auto norm = std::sqrt(normsq);
                        significant = (norm >= screening_threshold_);
                    }

                    if (significant) {
                        mx.lock();
                        ShellPairBlock block{{s1, s2}};
                        blocks.push_back(block);
                        mx.unlock();
                    }
                }
            }
        }  // end of compute

        // resort shell list in increasing order
        std::sort(blocks.begin(), blocks.end(), [](auto& l, auto& r) {
            return l[0].first == r[0].first ? l[0].second < r[0].second : l[0].first < r[0].first;
        });

        // compute shellpair data assuming that we are computing to default_epsilon
        // N.B. only parallelized over 1 shell index
        size_t npairs = blocks.size();
        ShellPairData spdata(blocks.size());
#pragma omp parallel
        {
            int thread_id = 0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif
            for (int pair = 0; pair < npairs; ++pair) {
                if (pair % nthreads == thread_id) {
                    auto s1 = blocks[pair][0].first;
                    auto s2 = blocks[pair][0].second;
                    spdata[pair] = std::make_shared<libint2::ShellPair>(bs1->l2_shell(s1), bs2->l2_shell(s2),
                                                                        std::log(screening_threshold_));
                }
            }
        }

        return std::make_tuple(blocks, spdata);
    }

    //! Shell pair information
    ShellPairData pairs12_, pairs34_;

    //! The threshold below which shell pairs are neglected on the basis of their overlap
    double screening_threshold_;

   public:
    //! Constructor. Use an IntegralFactory to create this object.
    Libint2TwoElectronInt(const IntegralFactory* integral, int deriv = 0, double screening_threshold = 0,
                          bool use_shell_pairs = false)
        : TwoBodyAOInt(integral, deriv), screening_threshold_(screening_threshold) {
        // Initialize libint static data
        libint2::initialize();

        // Figure out some information to initialize libint/libderiv with
        // 1. Maximum angular momentum
        int max_am = std::max(std::max(basis1()->max_am(), basis2()->max_am()),
                              std::max(basis3()->max_am(), basis4()->max_am()));
        // 2. Maximum number of primitive combinations
        int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                                 std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
        // TODO figure out how this threshold should be set, depending on the coefficients.  There are some examples in
        // Ed's tests
        libint2_ = libint2::Engine(op, max_nprim, max_am, deriv);
        libint2_.set(bk);

        int num_chunks;
        switch (deriv) {
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

        results_.resize(num_chunks * basis1()->max_function_per_shell() * basis2()->max_function_per_shell() *
                        basis3()->max_function_per_shell() * basis4()->max_function_per_shell());
        target_full_ = results_.data();
        target_ = target_full_;

        // These kinds of checks should be implemented when we figure out how to query the libint2 API
        // Make sure libint is compiled to handle our max AM
        // if (max_am >= LIBINT_MAX_AM) {
        //    outfile->Printf(
        //        "ERROR: ERI - Libint cannot handle angular momentum this high (%d).\n"
        //        "       Rebuild Libint with MAX_AM_ERI at least %d.\n",
        //        max_am, max_am);
        //    throw LimitExceeded<int>(
        //        "ERI - Libint cannot handle angular momentum this high.\n"
        //        "Rebuild Libint with MAX_AM_ERI at least (actual).\n",
        //        LIBINT_MAX_AM - 1, max_am, __FILE__, __LINE__);
        //} else if (deriv_ == 1 && max_am >= LIBDERIV_MAX_AM1) {
        //    outfile->Printf(
        //        "ERROR: ERI - Libint cannot handle angular momentum this high (%d) for first derivatives.\n"
        //        "     Rebuild Libint with MAX_AM_ERI at least %d.\n",
        //        max_am, max_am + 1);
        //    throw LimitExceeded<int>(
        //        "ERI - Libint cannot handle angular momentum this high.\n"
        //        "Rebuild Libint with MAX_AM_ERI at least (actual + 1).\n",
        //        LIBDERIV_MAX_AM1 - 1, max_am, __FILE__, __LINE__);
        //} else if (deriv_ == 2 && max_am >= LIBDERIV_MAX_AM12) {
        //    outfile->Printf(
        //        "ERROR: ERI - Libint cannot handle angular momentum this high (%d) for second derivatives.\n"
        //        "       Reconfigure Libint with MAX_AM_ERI at least %d\n",
        //        max_am, max_am + 2);
        //    throw LimitExceeded<int>(
        //        "ERI - Libint cannot handle angular momentum this high.\n"
        //        "Rebuild Libint with MAX_AM_ERI at least (actual + 2).\n",
        //        LIBDERIV_MAX_AM12 - 1, max_am, __FILE__, __LINE__);
        //} else if (deriv_ > 2) {
        //    outfile->Printf("ERROR: ERI - Cannot compute higher than second derivatives.");
        //    throw PSIEXCEPTION("ERI - Cannot compute higher than second derivatives.");
        //}

        // We want to allow shell pair usage, even if not all basis sets are the same.  For example
        // the (mn|P) integrals in density fitting should still use the shell pairs in the bra.
        if (use_shell_pairs) {
            std::tie(blocks12_, pairs12_) = build_shell_pair_list(basis1(), basis2());
            if (basis1() == basis3() && basis2() == basis4()) {
                pairs34_ = pairs12_;
                blocks34_ = blocks34_;
            } else {
                std::tie(blocks34_, pairs34_) = build_shell_pair_list(basis3(), basis4());
            }
        }
    }

    ~Libint2TwoElectronInt() override { libint2::finalize(); }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator& shellIter) override {
        return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
    }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(int s1, int s2, int s3, int s4) override {
#ifdef MINTS_TIMER
        timer_on("Libint2ERI::compute_shell");
#endif

        if (force_cartesian_) {
            throw PSIEXCEPTION("TwoElectronInt: bad instruction routing.");
        }
        int n1 = bs1_->shell(s1).nfunction();
        int n2 = bs2_->shell(s2).nfunction();
        int n3 = bs3_->shell(s3).nfunction();
        int n4 = bs4_->shell(s4).nfunction();
        size_t ntot = n1 * n2 * n3 * n4;

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        libint2_.compute2<op, bk, 0l>(sh1, sh2, sh3, sh4);

        const auto res = libint2_.results()[0];
        if (res == nullptr) {
            std::fill_n(target_, ntot, 0.0);
            return 0;
        }

        std::copy_n(res, ntot, target_);

#ifdef MINTS_TIMER
        timer_off("Libint2ERI::compute_shell");
#endif
        return ntot;
    }

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv1(int s1, int s2, int s3, int s4) override {
#ifdef MINTS_TIMER
        timer_on("Libint2ERI::compute_shell_deriv1");
#endif

        if (force_cartesian_) {
            throw PSIEXCEPTION("TwoElectronInt: bad instruction routing.");
        }
        int n1 = bs1_->shell(s1).nfunction();
        int n2 = bs2_->shell(s2).nfunction();
        int n3 = bs3_->shell(s3).nfunction();
        int n4 = bs4_->shell(s4).nfunction();
        size_t ntot = 12 * n1 * n2 * n3 * n4;
        printf("Ntot is %d    %d %d %d %d\n", ntot, n1, n2, n3, n4);

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        libint2_.compute2<op, bk, 1l>(sh1, sh2, sh3, sh4);

        const auto res = libint2_.results()[0];
        if (res == nullptr) {
            std::fill_n(target_, ntot, 0.0);
            return 0;
        }

        for (int i = 0; i < ntot; ++i) printf("%10.6f ", res[i]);
        printf("\n");
        std::copy_n(res, ntot, target_);

#ifdef MINTS_TIMER
        timer_off("Libint2ERI::compute_shell_deriv1");
#endif
        return ntot;
    }

    /// Compute ERI second derivatives between 4 sheels. Result is stored in buffer.
    size_t compute_shell_deriv2(int, int, int, int) override {
        throw PSIEXCEPTION("Libint2 deriv2 NYI");
        return 1;
    }

    void compute_shell_blocks(int shellpair12, int shellpair34, int npair12 = -1, int npair34 = -1) override {
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

        const auto* sp12 = pairs12_[shellpair12].get();
        const auto* sp34 = pairs34_[shellpair34].get();

        int n1 = bs1_->shell(s1).nfunction();
        int n2 = bs2_->shell(s2).nfunction();
        int n3 = bs3_->shell(s3).nfunction();
        int n4 = bs4_->shell(s4).nfunction();

        size_t ntot = n1 * n2 * n3 * n4;

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        libint2_.compute2<op, bk, 0L>(sh1, sh2, sh3, sh4, sp12, sp34);

        const auto res = libint2_.results()[0];
        if (res == nullptr) {
            std::fill_n(target_, ntot, 0.0);
        } else {
            std::copy_n(res, ntot, target_);
        }

#ifdef MINTS_TIMER
        timer_off("Libint2ERI::compute_shell_blocks");
#endif
    }
};

class ERI : public TwoElectronInt {
   public:
    ERI(const IntegralFactory* integral, int deriv = 0, bool use_shell_pairs = false);
    ~ERI() override;
};

class F12 : public TwoElectronInt {
   public:
    F12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv = 0,
        bool use_shell_pairs = false);
    ~F12() override;
};

class F12Scaled : public TwoElectronInt {
   public:
    F12Scaled(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv = 0,
              bool use_shell_pairs = false);
    ~F12Scaled() override;
};

class F12Squared : public TwoElectronInt {
   public:
    F12Squared(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv = 0,
               bool use_shell_pairs = false);
    ~F12Squared() override;
};

class F12G12 : public TwoElectronInt {
   public:
    F12G12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv = 0,
           bool use_shell_pairs = false);
    ~F12G12() override;
};

class F12DoubleCommutator : public TwoElectronInt {
   public:
    F12DoubleCommutator(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv = 0,
                        bool use_shell_pairs = false);
    ~F12DoubleCommutator() override;
};

class ErfERI : public TwoElectronInt {
   public:
    ErfERI(double omega, const IntegralFactory* integral, int deriv = 0, bool use_shell_pairs = false);
    ~ErfERI() override;

    void setOmega(double omega);
};

class ErfComplementERI : public TwoElectronInt {
   public:
    ErfComplementERI(double omega, const IntegralFactory* integral, int deriv = 0, bool use_shell_pairs = false);
    ~ErfComplementERI() override;

    void setOmega(double omega);
};

}  // namespace psi

#endif
