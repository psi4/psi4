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

#ifndef LIBMINTS_INTEGRALPARAMETERS_H
#define LIBMINTS_INTEGRALPARAMETERS_H

namespace psi {

class Vector;

class PSI_API IntegralParameters {
   private:
    size_t nparam_;

   public:
    IntegralParameters(size_t nparam = 0) : nparam_(nparam) {}
    virtual ~IntegralParameters() {}

    size_t nparam() const { return nparam_; }
};

class PSI_API CorrelationFactor : public IntegralParameters {
   private:
    double *coeff_;
    double *exponent_;

   public:
    CorrelationFactor(size_t nparam);
    CorrelationFactor(std::shared_ptr<Vector> coeff, std::shared_ptr<Vector> exponent);
    ~CorrelationFactor() override;

    virtual double slater_exponent() const { return 1.0; }
    void set_params(std::shared_ptr<Vector> coeff, std::shared_ptr<Vector> exponent);
    double *exponent() const { return exponent_; }
    double *coeff() const { return coeff_; }
};

class PSI_API FittedSlaterCorrelationFactor : public CorrelationFactor {
   private:
    double slater_exponent_;

   public:
    double slater_exponent() const override { return slater_exponent_; }

    FittedSlaterCorrelationFactor(double exponent);
    double exponent() { return slater_exponent_; }
};

}  // namespace psi

#endif  // INTEGRALPARAMETERS_H
