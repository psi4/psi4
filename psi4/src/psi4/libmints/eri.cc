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

#include "psi4/libmints/eri.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/fjt.h"
#include "psi4/libqt/qt.h"

#include <libint2/shell.h>
#include <libint2/engine.h>
using namespace psi;

/////////
// Normal two-electron repulsion integrals
/////////

#ifdef ENABLE_Libint1
ERI::ERI(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    // The +1 is needed for derivatives to work.
    fjt_ = new Taylor_Fjt(
        basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1, 1e-15);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

ERI::~ERI() { delete fjt_; }

/////////
// F12
/////////

F12::F12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    fjt_ = new F12Fundamental(
        cf, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

F12::~F12() { delete fjt_; }

/////////
// F12Scaled
/////////

F12Scaled::F12Scaled(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv,
                     bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    fjt_ = new F12ScaledFundamental(
        cf, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

F12Scaled::~F12Scaled() { delete fjt_; }

/////////
// F12 squared
/////////

F12Squared::F12Squared(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv,
                       bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    fjt_ = new F12SquaredFundamental(
        cf, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

F12Squared::~F12Squared() { delete fjt_; }

/////////
// F12G12
/////////

F12G12::F12G12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    fjt_ = new F12G12Fundamental(
        cf, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

F12G12::~F12G12() { delete fjt_; }

/////////
// F12DoubleCommutator
/////////

F12DoubleCommutator::F12DoubleCommutator(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral,
                                         int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    fjt_ = new F12DoubleCommutatorFundamental(
        cf, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

F12DoubleCommutator::~F12DoubleCommutator() { delete fjt_; }

/////////
// ErfERI
/////////

ErfERI::ErfERI(double omega, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfFundamental(
        omega, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

ErfERI::~ErfERI() { delete fjt_; }

void ErfERI::setOmega(double omega) { (static_cast<ErfFundamental *>(fjt_))->setOmega(omega); }

/////////
// ErfComplementERI
/////////

ErfComplementERI::ErfComplementERI(double omega, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs) {
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfComplementFundamental(
        omega, basis1()->max_am() + basis2()->max_am() + basis3()->max_am() + basis4()->max_am() + deriv_ + 1);
    // form the blocking. We use the default
    setup_sieve();
    TwoBodyAOInt::create_blocks();
}

ErfComplementERI::~ErfComplementERI() { delete fjt_; }

void ErfComplementERI::setOmega(double omega) { (static_cast<ErfComplementFundamental *>(fjt_))->setOmega(omega); }
#endif  // ENABLE_Libint1

//// Libint2 implementation
Libint2ERI::Libint2ERI(const IntegralFactory *integral, double screening_threshold, int deriv, bool use_shell_pairs,
                       bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2ERI::Libint2ERI");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision =
        0.;  // equivalent in accuracy and timings is `std::numeric_limits<double>::epsilon() * 1e-30;`

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

    for (int der = 0; der <= deriv; ++der) {
        // set braket upon engine construction so particular LIBINT2_MAX_AM_eri[|2|3] limit governs validity
        engines_.emplace_back(libint2::Operator::coulomb, max_nprim, max_am, der, max_precision,
                              libint2::operator_traits<libint2::Operator::coulomb>::default_params(), braket_);
    }
    // set max_am for primary basis to be sieved, not all basis1234
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ =
        libint2::Engine(libint2::Operator::coulomb, max_nprim, max_am, 0, max_precision,
                        libint2::operator_traits<libint2::Operator::coulomb>::default_params(), libint2::BraKet::xx_xx);
    common_init();
    timer_off("Libint2ERI::Libint2ERI");
}

void Libint2ERI::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                  const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ERI::libint2wrapper0");
    }
}

void Libint2ERI::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                  const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ERI::libint2wrapper1");
    }
}

void Libint2ERI::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                  const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::coulomb, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ERI::libint2wrapper2");
    }
}

Libint2ERI::~Libint2ERI(){}

Libint2ErfERI::Libint2ErfERI(double omega, const IntegralFactory *integral, double screening_threshold, int deriv,
                             bool use_shell_pairs, bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2ErfERI::Libint2ErfERI");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision = 0.;

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

    for (int der = 0; der <= deriv; ++der) {
        engines_.emplace_back(libint2::Operator::erf_coulomb, max_nprim, max_am, der, max_precision, omega, braket_);
    }
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ = libint2::Engine(libint2::Operator::erf_coulomb, max_nprim, max_am, 0, max_precision, omega,
                                      libint2::BraKet::xx_xx);
    common_init();
    timer_off("Libint2ErfERI::Libint2ErfERI");
}

void Libint2ErfERI::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                     const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                     const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ErfERI::libint2wrapper0");
    }
}

void Libint2ErfERI::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                     const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                     const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ErfERI::libint2wrapper1");
    }
}
void Libint2ErfERI::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                     const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                     const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::erf_coulomb, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                            sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ErfERI::libint2wrapper2");
    }
}

Libint2ErfERI::~Libint2ErfERI(){}

Libint2ErfComplementERI::Libint2ErfComplementERI(double omega, const IntegralFactory *integral,
                                                 double screening_threshold, int deriv, bool use_shell_pairs,
                                                 bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2ErfComplementERI::Libint2ErfComplementERI");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision = 0.;

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

    for (int der = 0; der <= deriv; ++der) {
        engines_.emplace_back(libint2::Operator::erfc_coulomb, max_nprim, max_am, der, max_precision, omega, braket_);
    }
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ = libint2::Engine(libint2::Operator::erfc_coulomb, max_nprim, max_am, 0, max_precision, omega,
                                      libint2::BraKet::xx_xx);
    common_init();
    timer_off("Libint2ErfComplementERI::Libint2ErfComplementERI");
}

void Libint2ErfComplementERI::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2,
                                               const libint2::Shell &sh3, const libint2::Shell &sh4,
                                               const libint2::ShellPair *sp12, const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ErfComplementERI::libint2wrapper0");
    }
}

void Libint2ErfComplementERI::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2,
                                               const libint2::Shell &sh3, const libint2::Shell &sh4,
                                               const libint2::ShellPair *sp12, const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ErfComplementERI::libint2wrapper1");
    }
}
void Libint2ErfComplementERI::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2,
                                               const libint2::Shell &sh3, const libint2::Shell &sh4,
                                               const libint2::ShellPair *sp12, const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::erfc_coulomb, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                             sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2ErfComplementERI::libint2wrapper2");
    }
}

Libint2ErfComplementERI::~Libint2ErfComplementERI(){}

//// Libint2 implementation
Libint2YukawaERI::Libint2YukawaERI(double zeta, const IntegralFactory *integral, double screening_threshold, int deriv,
                                   bool use_shell_pairs, bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2YukawaERI::Libint2YukawaERI");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision =
        0.;  // equivalent in accuracy and timings is `std::numeric_limits<double>::epsilon() * 1e-30;`

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
        throw PSIEXCEPTION("Bad BraKet type in Libint2YukawaERI");
    }

    for (int der = 0; der <= deriv; ++der) {
        // set braket upon engine construction so particular LIBINT2_MAX_AM_eri[|2|3] limit governs validity
        engines_.emplace_back(libint2::Operator::yukawa, max_nprim, max_am, der, max_precision, zeta, braket_);
    }
    // set max_am for primary basis to be sieved, not all basis1234
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ =
        libint2::Engine(libint2::Operator::yukawa, max_nprim, max_am, 0, max_precision, zeta, libint2::BraKet::xx_xx);
    common_init();
    timer_on("Libint2YukawaERI::Libint2YukawaERI");
}

void Libint2YukawaERI::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                        const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                        const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::yukawa, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::yukawa, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::yukawa, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::yukawa, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2YukawaERI::libint2_wrapper0");
    }
}

void Libint2YukawaERI::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                        const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                        const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::yukawa, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::yukawa, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::yukawa, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::yukawa, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2YukawaERI::libint2_wrapper1");
    }
}
void Libint2YukawaERI::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                        const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                        const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::yukawa, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::yukawa, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::yukawa, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::yukawa, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2YukawaERI::libint2_wrapper2");
    }
}

Libint2YukawaERI::~Libint2YukawaERI(){}

/// F12

Libint2F12::Libint2F12(std::vector<std::pair<double, double>> exp_coeff, const IntegralFactory *integral,
                       double screening_threshold, int deriv, bool use_shell_pairs, bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2F12::Libint2F12");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision = 0.;

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

    for (int der = 0; der <= deriv; ++der) {
        engines_.emplace_back(libint2::Operator::cgtg, max_nprim, max_am, der, max_precision, exp_coeff, braket_);
    }
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ = libint2::Engine(libint2::Operator::cgtg, max_nprim, max_am, 0, max_precision, exp_coeff,
                                      libint2::BraKet::xx_xx);
    common_init();
    timer_off("Libint2F12::Libint2F12");
}

void Libint2F12::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                  const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::cgtg, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::cgtg, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::cgtg, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::cgtg, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12::libint2wrapper0");
    }
}

void Libint2F12::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                  const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::cgtg, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::cgtg, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::cgtg, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::cgtg, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12::libint2wrapper1");
    }
}
void Libint2F12::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                  const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                  const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::cgtg, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::cgtg, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::cgtg, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::cgtg, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12, sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12::libint2wrapper2");
    }
}

Libint2F12::~Libint2F12(){}

/// F12G12

Libint2F12G12::Libint2F12G12(std::vector<std::pair<double, double>> exp_coeff, const IntegralFactory *integral,
                             double screening_threshold, int deriv, bool use_shell_pairs, bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2F12G12::Libint2F12G12");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision = 0.;

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

    for (int der = 0; der <= deriv; ++der) {
        engines_.emplace_back(libint2::Operator::cgtg_x_coulomb, max_nprim, max_am, der, max_precision, exp_coeff,
                              braket_);
    }
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ = libint2::Engine(libint2::Operator::cgtg_x_coulomb, max_nprim, max_am, 0, max_precision, exp_coeff,
                                      libint2::BraKet::xx_xx);
    common_init();
    timer_off("Libint2F12G12::Libint2F12G12");
}

void Libint2F12G12::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                     const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                     const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12G12::libint2wrapper0");
    }
}

void Libint2F12G12::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                     const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                     const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12G12::libint2wrapper1");
    }
}
void Libint2F12G12::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2, const libint2::Shell &sh3,
                                     const libint2::Shell &sh4, const libint2::ShellPair *sp12,
                                     const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::cgtg_x_coulomb, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                               sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12G12::libint2wrapper2");
    }
}

Libint2F12G12::~Libint2F12G12(){}

/// F12DoubleCommutator

Libint2F12DoubleCommutator::Libint2F12DoubleCommutator(std::vector<std::pair<double, double>> exp_coeff,
                                                       const IntegralFactory *integral, double screening_threshold,
                                                       int deriv, bool use_shell_pairs, bool needs_exchange)
    : Libint2TwoElectronInt(integral, deriv, screening_threshold, use_shell_pairs, needs_exchange) {
    timer_on("Libint2F12DoubleCommutator::Libint2F12DoubleCommutator");
    int max_am =
        std::max(std::max(basis1()->max_am(), basis2()->max_am()), std::max(basis3()->max_am(), basis4()->max_am()));
    int max_nprim = std::max(std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive()),
                             std::max(basis3()->max_nprimitive(), basis4()->max_nprimitive()));
    const auto max_precision = 0.;

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

    for (int der = 0; der <= deriv; ++der) {
        engines_.emplace_back(libint2::Operator::delcgtg2, max_nprim, max_am, der, max_precision, exp_coeff, braket_);
    }
    max_am = bra_same_ ? basis1()->max_am() : ket_same_ ? basis3()->max_am() : 0;
    schwarz_engine_ = libint2::Engine(libint2::Operator::delcgtg2, max_nprim, max_am, 0, max_precision, exp_coeff,
                                      libint2::BraKet::xx_xx);
    common_init();
    timer_off("Libint2F12DoubleCommutator::Libint2F12DoubleCommutator");
}

void Libint2F12DoubleCommutator::libint2_wrapper0(const libint2::Shell &sh1, const libint2::Shell &sh2,
                                                  const libint2::Shell &sh3, const libint2::Shell &sh4,
                                                  const libint2::ShellPair *sp12, const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[0].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xx_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[0].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xs_xx, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[0].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xx_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[0].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xs_xs, 0>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12DoubleCommutator::libint2wrapper0");
    }
}

void Libint2F12DoubleCommutator::libint2_wrapper1(const libint2::Shell &sh1, const libint2::Shell &sh2,
                                                  const libint2::Shell &sh3, const libint2::Shell &sh4,
                                                  const libint2::ShellPair *sp12, const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[1].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xx_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[1].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xs_xx, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[1].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xx_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[1].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xs_xs, 1>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12DoubleCommutator::libint2wrapper1");
    }
}
void Libint2F12DoubleCommutator::libint2_wrapper2(const libint2::Shell &sh1, const libint2::Shell &sh2,
                                                  const libint2::Shell &sh3, const libint2::Shell &sh4,
                                                  const libint2::ShellPair *sp12, const libint2::ShellPair *sp34) {
    switch (braket_) {
        case libint2::BraKet::xx_xx:
            engines_[2].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xx_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xs_xx:
            engines_[2].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xs_xx, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xx_xs:
            engines_[2].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xx_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        case libint2::BraKet::xs_xs:
            engines_[2].compute2<libint2::Operator::delcgtg2, libint2::BraKet::xs_xs, 2>(sh1, sh2, sh3, sh4, sp12,
                                                                                         sp34);
            break;
        default:
            throw PSIEXCEPTION("Bad BraKet type in Libint2F12DoubleCommutator::libint2wrapper2");
    }
}

Libint2F12DoubleCommutator::~Libint2F12DoubleCommutator(){}
