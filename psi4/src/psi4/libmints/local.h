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

#ifndef _psi_src_lib_libmints_local_h_
#define _psi_src_lib_libmints_local_h_

#include <vector>
#include <memory>

namespace psi {
class Options;

class Matrix;

class BasisSet;

class Localizer
{
protected:

    // => Parameters <= //

    /// Print flag
    int print_;
    /// Debug flug
    int debug_;
    /// Bench flag
    int bench_;

    /// Relative convergence criteria
    double convergence_;
    /// Maximum number of iterations
    int maxiter_;

    /// Primary orbital basis set
    std::shared_ptr <BasisSet> primary_;
    /// Delocalized Orbitals
    std::shared_ptr <Matrix> C_;

    // => Targets <= //

    /// Localized Orbitals
    std::shared_ptr <Matrix> L_;
    /// MO -> LO transformation
    std::shared_ptr <Matrix> U_;
    /// Did the algorithm converge?
    bool converged_;

    /// Set defaults
    void common_init();

public:

    // => Constructors <= //

    Localizer(std::shared_ptr <BasisSet> primary_, std::shared_ptr <Matrix> C);

    virtual ~Localizer();

    static std::shared_ptr <Localizer> build(const std::string &type, std::shared_ptr <BasisSet> primary, std::shared_ptr <Matrix> C);

    static std::shared_ptr <Localizer> build(const std::string &type, std::shared_ptr <BasisSet> primary, std::shared_ptr <Matrix> C, Options &options);

    static std::shared_ptr <Localizer> build(std::shared_ptr <BasisSet> primary, std::shared_ptr <Matrix> C, Options &options);

    // => Computers <= //

    /// Print out the localization algorithm and parameters
    virtual void print_header() const = 0;

    /// Perform the localization algorithm
    virtual void localize() = 0;

    /// Given a Fock matrix in the original basis (usually diagonal), produce an ordered copy in the local basis, and reorder L and U
    std::shared_ptr <Matrix> fock_update(std::shared_ptr <Matrix> F_orig);

    // => Accessors <= //

    std::shared_ptr <Matrix> L() const
    { return L_; }

    std::shared_ptr <Matrix> U() const
    { return U_; }

    bool converged() const
    { return converged_; }

    // => Knobs <= //

    void set_print(int print)
    { print_ = print; }

    void set_debug(int debug)
    { debug_ = debug; }

    void set_bench(int bench)
    { bench_ = bench; }

    void set_convergence(double convergence)
    { convergence_ = convergence; }

    void set_maxiter(int maxiter)
    { maxiter_ = maxiter; }

};

class BoysLocalizer : public Localizer
{

protected:

    /// Set defaults
    void common_init();

public:
    BoysLocalizer(std::shared_ptr <BasisSet> primary, std::shared_ptr <Matrix> C);

    virtual ~BoysLocalizer();

    virtual void print_header() const;

    virtual void localize();

};

class PMLocalizer : public Localizer
{

protected:

    /// Set defaults
    void common_init();

public:
    PMLocalizer(std::shared_ptr <BasisSet> primary, std::shared_ptr <Matrix> C);

    virtual ~PMLocalizer();

    virtual void print_header() const;

    virtual void localize();

};

} //Namespace psi

#endif
