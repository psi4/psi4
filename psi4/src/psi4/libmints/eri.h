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

#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

#ifdef _OPENMP
#include <omp.h>
#endif
#include <numeric>
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/shellpair.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/psi4-dec.h"

namespace libint2 {
enum class BraKet;
class Engine;
class shell;
}  // namespace libint2
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

/// Libint2 implementation

/*! \ingroup MINTS
 *  \class Libint2TwoElectronInt
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class Libint2TwoElectronInt : public TwoBodyAOInt {
   protected:
    //! Libint2 engine
    std::vector<libint2::Engine> engines_;

    /// The engine to used to set up the Schwarz screening
    libint2::Engine schwarz_engine_;

    //! Shell pair information
    ShellPairData pairs12_, pairs34_;

    /// The type of shell combo to be handled by this object
    libint2::BraKet braket_;

    /// A vector of zeros that we can point to if libint2 gives us back a nullptr
    std::vector<double> zero_vec_;
    bool use_shell_pairs_;

    //! Setup metadata and screening info
    void common_init();

   public:
    //! Constructor. Use an IntegralFactory to create this object.
    Libint2TwoElectronInt(const IntegralFactory *integral, int deriv = 0, double screening_threshold = 0,
                          bool use_shell_pairs = false, bool needs_exchange = false);

    Libint2TwoElectronInt(const Libint2TwoElectronInt &rhs);

    ~Libint2TwoElectronInt() override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(const AOShellCombinationsIterator &shellIter) override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3, int s4,
                                   bool is_bra) override;

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    size_t compute_shell(int s1, int s2, int s3, int s4) override;

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv1(int s1, int s2, int s3, int s4) override;

    /// Compute ERI second derivatives between 4 shells. Result is stored in buffer.
    size_t compute_shell_deriv2(int s1, int s2, int s3, int s4) override;

    virtual void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                                  const libint2::ShellPair *sp34 = nullptr) = 0;
    virtual void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                                  const libint2::ShellPair *sp34 = nullptr) = 0;
    virtual void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                                  const libint2::ShellPair *sp34 = nullptr) = 0;

    void compute_shell_blocks(int shellpair12, int shellpair34, int npair12 = -1, int npair34 = -1) override;
};

class Libint2ERI : public Libint2TwoElectronInt {
   public:
    Libint2ERI(const IntegralFactory *integral, double screening_threshold, int deriv = 0, bool use_shell_pairs = false,
               bool needs_exchange = false);
    ~Libint2ERI() override;
    Libint2ERI *clone() const override { return new Libint2ERI(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

class Libint2ErfERI : public Libint2TwoElectronInt {
   public:
    Libint2ErfERI(double omega, const IntegralFactory *integral, double screening_threshold, int deriv = 0,
                  bool use_shell_pairs = false, bool needs_exchange = false);
    ~Libint2ErfERI() override;
    Libint2ErfERI *clone() const override { return new Libint2ErfERI(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

class Libint2ErfComplementERI : public Libint2TwoElectronInt {
   public:
    Libint2ErfComplementERI(double omega, const IntegralFactory *integral, double screening_threshold, int deriv = 0,
                            bool use_shell_pairs = false, bool needs_exchange = false);
    ~Libint2ErfComplementERI() override;
    Libint2ErfComplementERI *clone() const override { return new Libint2ErfComplementERI(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

class Libint2F12 : public Libint2TwoElectronInt {
   public:
    Libint2F12(std::vector<std::pair<double, double>> exp_coeff, const IntegralFactory *integral,
               double screening_threshold, int deriv = 0, bool use_shell_pairs = false, bool needs_exchange = false);
    ~Libint2F12() override;
    Libint2F12 *clone() const override { return new Libint2F12(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

inline std::vector<std::pair<double, double>> take_square(std::vector<std::pair<double, double>> input) {
    auto n = input.size();
    std::vector<std::pair<double, double>> output;
    for (int i = 0; i < n; ++i) {
        auto e_i = input[i].first;
        auto c_i = input[i].second;
        for (int j = i; j < n; ++j) {
            auto e_j = input[j].first;
            auto c_j = input[j].second;
            double scale = i == j ? 1.0 : 2.0;
            output.emplace_back(std::make_pair(e_i + e_j, scale * c_i * c_j));
        }
    }
    return output;
}

class Libint2F12Squared : public Libint2F12 {
   public:
    Libint2F12Squared(std::vector<std::pair<double, double>> exp_coeff, const IntegralFactory *integral,
                      double screening_threshold, int deriv = 0, bool use_shell_pairs = false,
                      bool needs_exchange = false)
        : Libint2F12(take_square(exp_coeff), integral, screening_threshold, deriv, use_shell_pairs, needs_exchange) {}
};

class Libint2F12G12 : public Libint2TwoElectronInt {
   public:
    Libint2F12G12(std::vector<std::pair<double, double>> exp_coeff, const IntegralFactory *integral,
                  double screening_threshold, int deriv = 0, bool use_shell_pairs = false, bool needs_exchange = false);
    ~Libint2F12G12() override;
    Libint2F12G12 *clone() const override { return new Libint2F12G12(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

class Libint2F12DoubleCommutator : public Libint2TwoElectronInt {
   public:
    Libint2F12DoubleCommutator(std::vector<std::pair<double, double>> exp_coeff, const IntegralFactory *integral,
                               double screening_threshold, int deriv = 0, bool use_shell_pairs = false,
                               bool needs_exchange = false);
    ~Libint2F12DoubleCommutator() override;
    Libint2F12DoubleCommutator *clone() const override { return new Libint2F12DoubleCommutator(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

class Libint2YukawaERI : public Libint2TwoElectronInt {
   public:
    Libint2YukawaERI(double zeta, const IntegralFactory *integral, double screening_threshold, int deriv = 0,
                     bool use_shell_pairs = false, bool needs_exchange = false);
    ~Libint2YukawaERI() override;
    Libint2YukawaERI *clone() const override { return new Libint2YukawaERI(*this); }

   protected:
    void libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
    void libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                          const libint2::Shell &sh4, const libint2::ShellPair *sp12 = nullptr,
                          const libint2::ShellPair *sp34 = nullptr) override;
};

}  // namespace psi

#endif
