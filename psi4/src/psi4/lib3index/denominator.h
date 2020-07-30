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

#ifndef three_index_denominator_H
#define three_index_denominator_H

#include "psi4/libmints/tensor.h"

namespace psi {

class Matrix;

// Denominator Factorizations (MP2-like for now)
class PSI_API Denominator {
   protected:
    // Denominator (w in rows, ia in column)
    SharedMatrix denominator_;

    // Pointer to active occupied orbital eigenvalues
    SharedVector_<double> eps_occ_;
    // Pointer to active virtual orbital eigenvalues
    SharedVector_<double> eps_vir_;
    // Number of vectors required to obtain given accuracy
    int nvector_;
    // Maximum error norm allowed in denominator
    double delta_;

    virtual void decompose() = 0;

   public:
    Denominator(SharedVector_<double> eps_occ, SharedVector_<double> eps_vir, double delta);
    virtual ~Denominator();

    // Factory method, algorithm should be LAPLACE or CHOLESKY
    static std::shared_ptr<Denominator> buildDenominator(const std::string& algorithm, SharedVector_<double> eps_occ,
                                                         SharedVector_<double> eps_vir, double delta);

    double delta() const { return delta_; }
    int nvector() const { return nvector_; }
    virtual void debug();
    SharedMatrix denominator() const { return denominator_; }
};

class PSI_API LaplaceDenominator : public Denominator {
   protected:
    // Fully split denominator (w in rows, i in columns)
    SharedMatrix denominator_occ_;
    // Fully split denominator (w in rows, a in columns)
    SharedMatrix denominator_vir_;

    void decompose() override;

   public:
    LaplaceDenominator(SharedVector_<double> eps_occ_, SharedVector_<double> eps_vir, double delta);
    ~LaplaceDenominator() override;
    void debug() override;
    SharedMatrix denominator_occ() const { return denominator_occ_; }
    SharedMatrix denominator_vir() const { return denominator_vir_; }
};

class PSI_API CholeskyDenominator : public Denominator {
   protected:
    void decompose() override;

   public:
    CholeskyDenominator(SharedVector_<double> eps_occ_, SharedVector_<double> eps_vir, double delta);
    ~CholeskyDenominator() override;
    void debug() override;
};

class PSI_API SAPTDenominator {
   protected:
    // Denominator (w in rows, ar in column) (monomer A)
    SharedMatrix denominatorA_;
    // Denominator (w in rows, bs in column) (monomer B)
    SharedMatrix denominatorB_;

    // Pointer to active occupied orbital eigenvalues (monomer A)
    SharedVector_<double> eps_occA_;
    // Pointer to active virtual orbital eigenvalues (monomer A)
    SharedVector_<double> eps_virA_;
    // Pointer to active occupied orbital eigenvalues (monomer B)
    SharedVector_<double> eps_occB_;
    // Pointer to active virtual orbital eigenvalues (monomer B)
    SharedVector_<double> eps_virB_;
    // Number of vectors required to obtain given accuracy
    int nvector_;
    // Maximum error norm allowed in denominator
    double delta_;
    // Crap all over the output file?
    bool debug_;

    virtual void decompose() = 0;
    void check_denom(SharedVector_<double>, SharedVector_<double>, SharedMatrix);

   public:
    SAPTDenominator(SharedVector_<double>, SharedVector_<double>, SharedVector_<double>, SharedVector_<double>, double,
                    bool);
    virtual ~SAPTDenominator();

    // Factory method, algorithm should be LAPLACE or CHOLESKY
    static std::shared_ptr<SAPTDenominator> buildDenominator(
        const std::string& algorithm, SharedVector_<double> eps_occA, SharedVector_<double> eps_virA,
        SharedVector_<double> eps_occB, SharedVector_<double> eps_virB, double delta, bool debug = false);

    double delta() const { return delta_; }
    int nvector() const { return nvector_; }
    virtual void debug();
    SharedMatrix denominatorA() const { return denominatorA_; }
    SharedMatrix denominatorB() const { return denominatorB_; }
};

class PSI_API SAPTLaplaceDenominator : public SAPTDenominator {
   protected:
    // Fully split denominator (w in rows, a in columns) (monomer A)
    SharedMatrix denominator_occA_;
    // Fully split denominator (w in rows, r in columns) (monomer A)
    SharedMatrix denominator_virA_;
    // Fully split denominator (w in rows, b in columns) (monomer B)
    SharedMatrix denominator_occB_;
    // Fully split denominator (w in rows, s in columns) (monomer B)
    SharedMatrix denominator_virB_;

    void decompose() override;
    void check_split(SharedVector_<double>, SharedVector_<double>, SharedMatrix, SharedMatrix);

   public:
    SAPTLaplaceDenominator(SharedVector_<double>, SharedVector_<double>, SharedVector_<double>, SharedVector_<double>,
                           double, bool debug = false);
    ~SAPTLaplaceDenominator() override;

    void debug() override;
    SharedMatrix denominator_occA() const { return denominator_occA_; }
    SharedMatrix denominator_virA() const { return denominator_virA_; }
    SharedMatrix denominator_occB() const { return denominator_occB_; }
    SharedMatrix denominator_virB() const { return denominator_virB_; }
};

class PSI_API SAPTCholeskyDenominator : public SAPTDenominator {
   protected:
    void decompose() override;

   public:
    SAPTCholeskyDenominator(SharedVector_<double>, SharedVector_<double>, SharedVector_<double>, SharedVector_<double>,
                            double, bool debug = false);
    ~SAPTCholeskyDenominator() override;
};

class PSI_API TLaplaceDenominator {
    // Pointer to active occupied orbital eigenvalues
    SharedVector_<double> eps_occ_;
    // Pointer to active virtual orbital eigenvalues
    SharedVector_<double> eps_vir_;
    // Maximum error norm allowed in denominator
    double delta_;

    // Fully split denominator (w in rows i in columns)
    SharedMatrix denominator_occ_;
    // Fully split denominator (w in rows i in columns)
    SharedMatrix denominator_vir_;
    // Number of vectors required to obtain given accuracy
    int nvector_;

    virtual void decompose();

   public:
    TLaplaceDenominator(SharedVector_<double> eps_occ, SharedVector_<double> eps_vir, double delta);
    virtual ~TLaplaceDenominator();

    double delta() const { return delta_; }
    int nvector() const { return nvector_; }
    virtual void debug();
    SharedMatrix denominator_occ() const { return denominator_occ_; }
    SharedMatrix denominator_vir() const { return denominator_vir_; }
};

}  // Namespace psi
#endif
