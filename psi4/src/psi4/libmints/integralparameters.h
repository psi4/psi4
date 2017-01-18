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

#ifndef LIBMINTS_INTEGRALPARAMETERS_H
#define LIBMINTS_INTEGRALPARAMETERS_H

namespace psi {

class Vector;

class IntegralParameters
{
private:
    unsigned int nparam_;
public:
    IntegralParameters(unsigned int nparam=0)
        : nparam_(nparam) { }
    virtual ~IntegralParameters() {}

    unsigned int nparam() const { return nparam_; }
};

class CorrelationFactor : public IntegralParameters
{
private:
    double *coeff_;
    double *exponent_;

public:
    CorrelationFactor(unsigned int nparam);
    CorrelationFactor(std::shared_ptr<Vector> coeff,
                      std::shared_ptr<Vector> exponent);
    virtual ~CorrelationFactor();

    virtual double slater_exponent() const { return 1.0; }
    void set_params(std::shared_ptr<Vector> coeff,
                    std::shared_ptr<Vector> exponent);
    double *exponent() const { return exponent_; }
    double *coeff() const { return coeff_; }
};

class FittedSlaterCorrelationFactor : public CorrelationFactor
{
private:
    double slater_exponent_;

public:

    virtual double slater_exponent() const { return slater_exponent_; }

    FittedSlaterCorrelationFactor(double exponent);
    double exponent(){return slater_exponent_;}
};

}

#endif // INTEGRALPARAMETERS_H
