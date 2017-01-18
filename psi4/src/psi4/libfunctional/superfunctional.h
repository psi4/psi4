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

#ifndef SUPERFUNCTIONAL_H
#define SUPERFUNCTIONAL_H

#include "psi4/libmints/typedefs.h"
#include <map>
#include <vector>
#include <cstdlib>
#include <string>
namespace psi {

class Options;
class Functional;

/**
 * SuperFunctional: High-level semilocal DFA object
 *
 * This class directly holds semilocal DFA objects
 * and sketches the definition of the generalized Kohn-Sham functional
 * via alpha and omega parameters for long-range and global hybrid
 * exact exchange mixing, and alpha and omega parameters for long-range
 * and global MP2-like correlation.
 *
 * A DFT functional is defined as:
 *
 * E_XC = E_X + E_C + E(-D)
 * E_X  = (1-\alpha_x) E_X^DFA [\omega_x] + \alpha_x E_X^HF + (1-\alpha_x) E_X^HF,LR
 * E_C  = E_C^DFA [\omega_c] + \alpha_c E_C^MP2
 *
 * If SCS-MP2:
 * E_C  = E_C^DFA [\omega_c] + \alpha_c_ss E_C^SS-MP2 + \alpha_c_os E_C^OS-MP2
 *
 * (Note: earlier versions of this code defined --- but did not fully implement --- a
 * range-separated MP2 correction, like \alpha_c E_C^MP2 + (1-\alpha_c) E_C^MP2,LR.
 * This is not currently supported given the reworking of how alpha_c works for
 * scaling the MP2 correlation (it needs to be independent of the DFT corrleation
 * in general, not related by \alpha_c and (1-\alpha_c).)
 *
 * E(-D) is an empirical dispersion correction of some form
 *
 **/
class SuperFunctional {

protected:

    // => Meta-Data <= //

    std::string name_;
    std::string description_;
    std::string citation_;

    // => Exchange-side DFA functionals <= //

    std::vector<std::shared_ptr<Functional> > x_functionals_;
    double x_alpha_;
    double x_omega_;

    // => Correlation-side DFA functionals <= //

    std::vector<std::shared_ptr<Functional> > c_functionals_;
    double c_alpha_;
    double c_ss_alpha_;
    double c_os_alpha_;
    double c_omega_;

    // => Functional values and partials <= //

    int max_points_;
    int deriv_;
    std::map<std::string, SharedVector> values_;

    // The omegas or alphas have changed, we're in a GKS environment.
    // Update the short-range DFAs
    void partition_gks();
    // Set up a null Superfunctional
    void common_init();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    SuperFunctional();
    virtual ~SuperFunctional();

    static std::shared_ptr<SuperFunctional> blank();

    // Allocate values (MUST be called after adding new functionals to the superfunctional)
    void allocate();

    // => Computers <= //

    std::map<std::string, SharedVector>& compute_functional(const std::map<std::string, SharedVector>& vals, int npoints = -1);
    void test_functional(SharedVector rho_a,
                         SharedVector rho_b,
                         SharedVector gamma_aa,
                         SharedVector gamma_ab,
                         SharedVector gamma_bb,
                         SharedVector tau_a,
                         SharedVector tau_b);

    // => Input/Output <= //

    std::map<std::string, SharedVector>& values() { return values_; }
    SharedVector value(const std::string& key);

    std::vector<std::shared_ptr<Functional> >& x_functionals() { return x_functionals_; }
    std::vector<std::shared_ptr<Functional> >& c_functionals() { return c_functionals_; }

    std::shared_ptr<Functional> x_functional(const std::string& name);
    std::shared_ptr<Functional> c_functional(const std::string& name);
    void add_x_functional(std::shared_ptr<Functional> fun);
    void add_c_functional(std::shared_ptr<Functional> fun);

    // => Setters <= //

    void set_name(const std::string & name) { name_ = name; }
    void set_description(const std::string & description) { description_ = description; }
    void set_citation(const std::string & citation) { citation_ = citation; }

    void set_max_points(int max_points) { max_points_ = max_points; allocate(); }
    void set_deriv(int deriv) { deriv_ = deriv;  allocate(); }

    void set_x_omega(double omega) { x_omega_ = omega; partition_gks(); }
    void set_c_omega(double omega) { c_omega_ = omega; partition_gks(); }
    void set_x_alpha(double alpha) { x_alpha_ = alpha; partition_gks(); }
    void set_c_alpha(double alpha) { c_alpha_ = alpha; partition_gks(); }
    void set_c_ss_alpha(double alpha) { c_ss_alpha_ = alpha; partition_gks(); }
    void set_c_os_alpha(double alpha) { c_os_alpha_ = alpha; partition_gks(); }

    // => Accessors <= //

    std::string name() const { return name_; }
    std::string description() const { return description_; }
    std::string citation() const { return citation_; }

    int ansatz() const;
    int max_points() const { return max_points_; }
    int deriv() const { return deriv_; }

    double x_omega() const { return x_omega_; }
    double c_omega() const { return c_omega_; }
    double x_alpha() const { return x_alpha_; }
    double c_alpha() const { return c_alpha_; }
    double c_ss_alpha() const { return c_ss_alpha_; }
    double c_os_alpha() const { return c_os_alpha_; }

    bool needs_xc() const { return ((c_functionals_.size() + x_functionals_.size()) > 0); }
    bool is_meta() const;
    bool is_gga() const;
    bool is_x_lrc() const { return x_omega_ != 0.0; }
    bool is_c_lrc() const { return c_omega_ != 0.0; }
    bool is_x_hybrid() const { return x_alpha_ != 0.0; }
    bool is_c_hybrid() const { return c_alpha_ != 0.0; }
    bool is_c_scs_hybrid() const { return c_os_alpha_ != 0.0 || c_ss_alpha_ != 0.0; }

    // => Utility <= //

    void print(std::string OutFileRMR = "outfile", int print = 1) const;
    void py_print() const { print("outfile", 1); }
    void py_print_detail(int level) const { print("outfile", level); }

};

}

#endif
