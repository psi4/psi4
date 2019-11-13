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

#ifdef _OPENMP
#include <omp.h>
#endif
#include <unordered_map>
#include <numeric>
#include <libint2/engine.h>
#include <libint2/shell.h>
#include <libint/libint.h>
#include <libderiv/libderiv.h>
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/shellpair.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/psi4-dec.h"

namespace psi {

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
struct L1ShellPair {
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

    //! The number of integrals in the current shell quartet
    int batchsize_;

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
    std::shared_ptr<std::vector<std::vector<L1ShellPair>>> pairs12_, pairs34_;

    //! Evaluates how much memory (in doubles) is needed to store shell pair data
    size_t memory_to_store_shell_pairs(const std::shared_ptr<BasisSet>&, const std::shared_ptr<BasisSet>&);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;
    size_t compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3, int s4, bool is_bra) override;

   public:
    //! Constructor. Use an IntegralFactory to create this object.
    TwoElectronInt(const IntegralFactory* integral, int deriv = 0, bool use_shell_pairs = false);

    ~TwoElectronInt() override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator&) override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(int s1, int s2, int s3, int s4) override;

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv1(int s1, int s2, int s3, int s4) override;

    /// Compute ERI second derivatives between 4 sheels. Result is stored in buffer.
    size_t compute_shell_deriv2(int s1, int s2, int s3, int s4) override;
};

/*! \ingroup MINTS
 *  \class Libint2TwoElectronInt
 *  \brief Capable of computing two-electron repulsion integrals.
 */
template <libint2::Operator op, libint2::BraKet bk>
class Libint2TwoElectronInt : public TwoBodyAOInt {
   protected:
    //! Libint2 engine
    std::vector<libint2::Engine> engines_;
    std::vector<double> results_;

    //! Shell pair information
    ShellPairData pairs12_, pairs34_;

    /// A vector of zeros that we can point to if libint2 gives us back a nullptr
    std::vector<double> zero_vec_;
    bool use_shell_pairs_;
    //! The threshold below which shell pairs are neglected on the basis of their overlap
    double screening_threshold_;

   public:
    //! Constructor. Use an IntegralFactory to create this object.
    Libint2TwoElectronInt(const IntegralFactory* integral, int deriv = 0, double screening_threshold = 0,
                          bool use_shell_pairs = false)
        : TwoBodyAOInt(integral, deriv), screening_threshold_(screening_threshold), use_shell_pairs_(use_shell_pairs) {
        // Initialize libint static data
        libint2::initialize();

        // Figure out some information to initialize libint/libderiv with
        // 1. Maximum angular momentum
        int max_am = std::max(std::max(basis1()->max_am(), basis2()->max_am()),
                              std::max(basis3()->max_am(), basis4()->max_am()));
        // 2. Maximum number of primitive combinations
        int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                                 std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
        for (int der = 0; der <= deriv; ++der) {
            engines_.emplace_back(op, max_nprim, max_am, der);
            engines_[der].set(bk);
        }

        zero_vec_ = std::vector<double>(basis1()->max_function_per_shell() * basis2()->max_function_per_shell() *
                                            basis3()->max_function_per_shell() * basis4()->max_function_per_shell(),
                                        0.0);
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
        buffers_.resize(num_chunks);

        target_full_ = const_cast<double*>(engines_[0].results()[0]);
        target_ = target_full_;

        // Make sure the engine can handle the type of integral used to build a sieve
        engines_[0].set(libint2::BraKet::xx_xx);
        setup_sieve();
        // Reset the engine type back to the general case needed
        engines_[0].set(bk);
        create_blocks();
        const auto max_engine_precision = std::numeric_limits<double>::epsilon() / 1e10;

        size_t npairs = shell_pairs_bra_.size();
        pairs12_.resize(npairs);
#pragma omp parallel for
        for (int pair = 0; pair < npairs; ++pair) {
            auto s1 = shell_pairs_bra_[pair].first;
            auto s2 = shell_pairs_bra_[pair].second;
            pairs12_[pair] = std::make_shared<libint2::ShellPair>(basis1()->l2_shell(s1), basis2()->l2_shell(s2),
                                                                  std::log(max_engine_precision));
        }
        npairs = shell_pairs_ket_.size();
        pairs34_.resize(npairs);
#pragma omp parallel for
        for (int pair = 0; pair < npairs; ++pair) {
            auto s3 = shell_pairs_ket_[pair].first;
            auto s4 = shell_pairs_ket_[pair].second;
            pairs34_[pair] = std::make_shared<libint2::ShellPair>(
                basis3()->l2_shell(s3), basis4()->l2_shell(s4), std::log(max_engine_precision));
        }
    }

    ~Libint2TwoElectronInt() override { libint2::finalize(); }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator& shellIter) override {
        return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
    }

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3, int s4, bool is_bra) override {
#ifdef MINTS_TIMER
        timer_on("Libint2ERI::compute_shell_for_sieve");
#endif

        auto sh1 = bs->l2_shell(s1);
        auto sh2 = bs->l2_shell(s2);
        auto sh3 = bs->l2_shell(s3);
        auto sh4 = bs->l2_shell(s4);

        engines_[0].compute2<op, libint2::BraKet::xx_xx, 0l>(sh1, sh2, sh3, sh4);

        size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();

        buffers_[0] = target_full_ = const_cast<double*>(engines_[0].results()[0]);
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

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(int s1, int s2, int s3, int s4) override {
#ifdef MINTS_TIMER
        timer_on("Libint2ERI::compute_shell");
#endif
        if (force_cartesian_) {
            throw PSIEXCEPTION("TwoElectronInt: bad instruction routing.");
        }

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        engines_[0].compute2<op, bk, 0l>(sh1, sh2, sh3, sh4);

        size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();

        buffers_[0] = target_full_ = const_cast<double*>(engines_[0].results()[0]);
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

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv1(int s1, int s2, int s3, int s4) override {
#ifdef MINTS_TIMER
        timer_on("Libint2ERI::compute_shell_deriv1");
#endif
        if (force_cartesian_) {
            throw PSIEXCEPTION("TwoElectronInt: bad instruction routing.");
        }

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        engines_[1].compute2<op, bk, 1l>(sh1, sh2, sh3, sh4);

        size_t shell_size = sh1.size() * sh2.size() * sh3.size() * sh4.size();
        size_t ntot = 0;

        for (int i = 0; i < 12; ++i) {
            if (engines_[1].results()[i]) {
                buffers_[i] = engines_[1].results()[i];
                ntot += shell_size;
            } else {
                buffers_[i] = zero_vec_.data();
            }
        }

#ifdef MINTS_TIMER
        timer_off("Libint2ERI::compute_shell_deriv1");
#endif
        return ntot;
    }

    /// Compute ERI second derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv2(int s1, int s2, int s3, int s4) override {
#ifdef MINTS_TIMER
        timer_on("Libint2ERI::compute_shell_deriv2");
#endif
        if (force_cartesian_) {
            throw PSIEXCEPTION("TwoElectronInt: bad instruction routing.");
        }

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        engines_[2].compute2<op, bk, 2l>(sh1, sh2, sh3, sh4);

        size_t shell_size = sh1.size() * sh2.size() * sh3.size() * sh4.size();
        size_t ntot = 0;

        for (int i = 0; i < 78; ++i) {
            if (engines_[2].results()[i]) {
                buffers_[i] = engines_[2].results()[i];
                ntot += shell_size;
            } else {
                buffers_[i] = zero_vec_.data();
            }
        }

#ifdef MINTS_TIMER
        timer_off("Libint2ERI::compute_shell_deriv2");
#endif
        return ntot;
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

        auto sh1 = bs1_->l2_shell(s1);
        auto sh2 = bs2_->l2_shell(s2);
        auto sh3 = bs3_->l2_shell(s3);
        auto sh4 = bs4_->l2_shell(s4);

        size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();

        //if(bra_same_ && ket_same_) {
            const auto* sp12 = pairs12_[shellpair12].get();
            const auto* sp34 = pairs34_[shellpair34].get();
            engines_[0].compute2<op, bk, 0L>(sh1, sh2, sh3, sh4, sp12, sp34);
        //} else {
        //    engines_[0].compute2<op, bk, 0L>(sh1, sh2, sh3, sh4);
        //}

        target_full_ = const_cast<double*>(engines_[0].results()[0]);
        if (target_full_) {
            buffers_[0] = engines_[0].results()[0];
        } else {
            target_full_ = zero_vec_.data();
            ntot = 0;
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
