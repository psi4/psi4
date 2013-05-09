/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#ifndef _psi_src_lib_libmints_local_h_
#define _psi_src_lib_libmints_local_h_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

namespace psi {

class Matrix;
class BasisSet;

class Localizer {

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
    boost::shared_ptr<BasisSet> primary_;
    /// Delocalized Orbitals
    boost::shared_ptr<Matrix> C_;
    
    // => Targets <= //
    
    /// Localized Orbitals 
    boost::shared_ptr<Matrix> L_;
    /// MO -> LO transformation 
    boost::shared_ptr<Matrix> U_;
    /// Did the algorithm converge?
    bool converged_;

    /// Set defaults
    void common_init();

public:

    // => Constructors <= //

    Localizer(boost::shared_ptr<BasisSet> primary_, boost::shared_ptr<Matrix> C);
    
    virtual ~Localizer();

    static boost::shared_ptr<Localizer> build(const std::string& type, boost::shared_ptr<BasisSet> primary, boost::shared_ptr<Matrix> C, Options& options);
    static boost::shared_ptr<Localizer> build(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<Matrix> C, Options& options);
    
    // => Computers <= //

    virtual void print_header() const = 0;
    virtual void localize() = 0;

    // => Accessors <= //

    boost::shared_ptr<Matrix> L() const { return L_; }
    boost::shared_ptr<Matrix> U() const { return U_; }
    bool converged() const { return converged_; }

    // => Knobs <= //

    void set_print(int print) { print_ = print; }
    void set_debug(int debug) { debug_ = debug; }
    void set_bench(int bench) { bench_ = bench; }
    void set_convergence(double convergence) { convergence_ = convergence; }
    void set_maxiter(int maxiter) { maxiter_ = maxiter; }

};

class BoysLocalizer : public Localizer {

protected:
    
    /// Set defaults
    void common_init();

public:
    BoysLocalizer(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<Matrix> C);

    virtual ~BoysLocalizer();

    virtual void print_header() const;
    virtual void localize();

};

class PMLocalizer : public Localizer {

protected:
    
    /// Set defaults
    void common_init();

public:
    PMLocalizer(boost::shared_ptr<BasisSet> primary, boost::shared_ptr<Matrix> C);

    virtual ~PMLocalizer();

    virtual void print_header() const;
    virtual void localize();

};

} //Namespace psi

#endif
