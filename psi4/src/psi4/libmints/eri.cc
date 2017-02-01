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

#include "psi4/libmints/eri.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/fjt.h"
;
using namespace psi;

/////////
// Normal two-electron repulsion integrals
/////////

ERI::ERI(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new Taylor_Fjt(basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1, 1e-15);
}

ERI::~ERI()
{
    delete fjt_;
}

/////////
// F12
/////////

F12::F12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    fjt_ = new F12Fundamental(cf,
                              basis1()->max_am() +
                              basis2()->max_am() +
                              basis3()->max_am() +
                              basis4()->max_am() +
                              deriv_+1);
}

F12::~F12()
{
    delete fjt_;
}

/////////
// F12Scaled
/////////

F12Scaled::F12Scaled(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    fjt_ = new F12ScaledFundamental(cf,
                              basis1()->max_am() +
                              basis2()->max_am() +
                              basis3()->max_am() +
                              basis4()->max_am() +
                              deriv_+1);
}

F12Scaled::~F12Scaled()
{
    delete fjt_;
}

/////////
// F12 squared
/////////

F12Squared::F12Squared(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    fjt_ = new F12SquaredFundamental(cf,
                                     basis1()->max_am() +
                                     basis2()->max_am() +
                                     basis3()->max_am() +
                                     basis4()->max_am() +
                                     deriv_+1);
}

F12Squared::~F12Squared()
{
    delete fjt_;
}

/////////
// F12G12
/////////

F12G12::F12G12(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    fjt_ = new F12G12Fundamental(cf,
                                 basis1()->max_am() +
                                 basis2()->max_am() +
                                 basis3()->max_am() +
                                 basis4()->max_am() +
                                 deriv_+1);
}

F12G12::~F12G12()
{
    delete fjt_;
}

/////////
// F12DoubleCommutator
/////////

F12DoubleCommutator::F12DoubleCommutator(std::shared_ptr<CorrelationFactor> cf, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    fjt_ = new F12DoubleCommutatorFundamental(cf,
                                              basis1()->max_am() +
                                              basis2()->max_am() +
                                              basis3()->max_am() +
                                              basis4()->max_am() +
                                              deriv_+1);
}

F12DoubleCommutator::~F12DoubleCommutator()
{
    delete fjt_;
}

/////////
// ErfERI
/////////

ErfERI::ErfERI(double omega, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfFundamental(omega,
                          basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1);
}

ErfERI::~ErfERI()
{
    delete fjt_;
}

void ErfERI::setOmega(double omega)
{
    (static_cast<ErfFundamental*>(fjt_))->setOmega(omega);
}

/////////
// ErfComplementERI
/////////

ErfComplementERI::ErfComplementERI(double omega, const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoElectronInt(integral, deriv, use_shell_pairs)
{
    // The +1 is needed for derivatives to work.
    fjt_ = new ErfComplementFundamental(omega,
                          basis1()->max_am() +
                          basis2()->max_am() +
                          basis3()->max_am() +
                          basis4()->max_am() +
                          deriv_+1);
}

ErfComplementERI::~ErfComplementERI()
{
    delete fjt_;
}

void ErfComplementERI::setOmega(double omega)
{
    (static_cast<ErfComplementFundamental*>(fjt_))->setOmega(omega);
}
