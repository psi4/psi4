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

#include "psi4/libmints/vector.h"
#include "psi4/libmints/integralparameters.h"

using namespace psi;

CorrelationFactor::CorrelationFactor(unsigned int nparam)
    : IntegralParameters(nparam)
{
}

CorrelationFactor::CorrelationFactor(std::shared_ptr<Vector> coeff, std::shared_ptr<Vector> exponent)
    : IntegralParameters(coeff->dim())
{
    set_params(coeff, exponent);
}

CorrelationFactor::~CorrelationFactor()
{
    delete[] coeff_;
    delete[] exponent_;
}

void CorrelationFactor::set_params(std::shared_ptr<Vector> coeff, std::shared_ptr<Vector> exponent)
{
    int nparam = coeff->dim();
    if (nparam) {
        coeff_ = new double[nparam];
        exponent_ = new double[nparam];
        for (int i=0; i<nparam; ++i) {
            coeff_[i] = coeff->get(0, i);
            exponent_[i] = exponent->get(0, i);
        }
    }
}

FittedSlaterCorrelationFactor::FittedSlaterCorrelationFactor(double exponent)
    : CorrelationFactor(6)
{
    // Perform the fit.
    SharedVector exps(new Vector(6));
    SharedVector coeffs(new Vector(6));

    slater_exponent_ = exponent;

    // The fitting coefficients
    coeffs->set(0, 0, -0.3144);
    coeffs->set(0, 1, -0.3037);
    coeffs->set(0, 2, -0.1681);
    coeffs->set(0, 3, -0.09811);
    coeffs->set(0, 4, -0.06024);
    coeffs->set(0, 5, -0.03726);

    // and the exponents
    exps->set(0, 0, 0.2209);
    exps->set(0, 1, 1.004);
    exps->set(0, 2, 3.622);
    exps->set(0, 3, 12.16);
    exps->set(0, 4, 45.87);
    exps->set(0, 5, 254.4);

    // They just need to be scaled
    double expsq = exponent * exponent;
    exps->scale(expsq);
    set_params(coeffs, exps);
}
