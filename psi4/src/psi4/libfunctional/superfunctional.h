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

#ifndef SUPERFUNCTIONAL_H
#define SUPERFUNCTIONAL_H

#include "psi4/libmints/typedefs.h"
#include <map>
#include <vector>
#include <cstdlib>
#include <string>
namespace psi {

class Functional;
class BlockOPoints;

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
    bool locked_;

    // => Exchange-side DFA functionals <= //
    std::vector<std::shared_ptr<Functional> > x_functionals_;
    double x_alpha_;
    double x_beta_;
    double x_omega_;

    // => Correlation-side DFA functionals <= //
    std::vector<std::shared_ptr<Functional> > c_functionals_;
    double c_alpha_;
    double c_ss_alpha_;
    double c_os_alpha_;
    double c_omega_;

    // => Asymptotic corrections <= //
    bool needs_grac_;
    std::shared_ptr<Functional> grac_x_functional_;
    std::shared_ptr<Functional> grac_c_functional_;
    double grac_shift_;
    double grac_alpha_;
    double grac_beta_;

    // => VV10 parameters corrections <= //
    bool needs_vv10_;
    double vv10_b_;
    double vv10_c_;
    double vv10_beta_;

    // => Functional values and partials <= //
    bool libxc_xc_func_;
    int max_points_;
    int deriv_;
    std::map<std::string, SharedVector> values_;
    std::map<std::string, SharedVector> ac_values_;
    std::map<std::string, SharedVector> vv_values_;

    // Set up a null Superfunctional
    void common_init();

    // Check if we can edit this Superfunctional
    void can_edit();

public:

    // => Constructors (Use the factory constructor, or really know what's up) <= //

    SuperFunctional();
    virtual ~SuperFunctional();

    // Build a blank superfunctional
    static std::shared_ptr<SuperFunctional> blank();
    static std::shared_ptr<SuperFunctional> XC_build(std::string name, bool unpolarized);

    // Builds a worker version of the superfunctional
    std::shared_ptr<SuperFunctional> build_worker();

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

    // Compute the cache data for VV10 dispersion
    std::map<std::string, SharedVector> compute_vv10_cache(
        const std::map<std::string, SharedVector>& vals, std::shared_ptr<BlockOPoints> block,
        double rho_thresh, int npoints = -1, bool internal = false);

    // Copmutes the Cache data for VV10 dispersion
    double compute_vv10_kernel(const std::map<std::string, SharedVector>& vals,
                               const std::vector<std::map<std::string, SharedVector>>& vv10_cache,
                               std::shared_ptr<BlockOPoints> block, int npoints = -1);

    // => Input/Output <= //

    std::map<std::string, SharedVector>& values() { return values_; }
    SharedVector value(const std::string& key);

    std::vector<std::shared_ptr<Functional> >& x_functionals() { return x_functionals_; }
    std::vector<std::shared_ptr<Functional> >& c_functionals() { return c_functionals_; }
    std::shared_ptr<Functional> grac_x_functional() { return grac_x_functional_; }
    std::shared_ptr<Functional> grac_c_functional() { return grac_c_functional_; }

    std::shared_ptr<Functional> x_functional(const std::string& name);
    std::shared_ptr<Functional> c_functional(const std::string& name);
    void add_x_functional(std::shared_ptr<Functional> fun);
    void add_c_functional(std::shared_ptr<Functional> fun);
    void set_grac_x_functional(std::shared_ptr<Functional> fun) {
        needs_grac_ = true;
        grac_x_functional_ = fun;
    }
    void set_grac_c_functional(std::shared_ptr<Functional> fun) {
        needs_grac_ = true;
        grac_c_functional_ = fun;
    }

    // => Setters <= //

    void set_lock(bool locked) { locked_ = locked; }
    void set_name(const std::string & name) { name_ = name; }
    void set_description(const std::string & description) { description_ = description; }
    void set_citation(const std::string & citation) { citation_ = citation; }

    void set_max_points(int max_points) { max_points_ = max_points; }
    void set_deriv(int deriv) { deriv_ = deriv; }

    void set_x_omega(double omega);
    void set_c_omega(double omega);
    void set_x_alpha(double alpha);
    void set_x_beta(double beta);
    void set_c_alpha(double alpha);
    void set_c_ss_alpha(double alpha);
    void set_c_os_alpha(double alpha);
    void set_vv10_b(double b);
    void set_vv10_c(double c);
    void set_grac_shift(double grac_shift);
    void set_grac_alpha(double grac_alpha);
    void set_grac_beta(double grac_beta);

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
    double x_beta() const { return x_beta_; }
    double c_alpha() const { return c_alpha_; }
    double c_ss_alpha() const { return c_ss_alpha_; }
    double c_os_alpha() const { return c_os_alpha_; }
    double vv10_b() const { return vv10_b_; }
    double vv10_c() const { return vv10_c_; }
    double grac_shift() const { return grac_shift_; }
    double grac_alpha() const { return grac_alpha_; }
    double grac_beta() const { return grac_beta_; }

    bool needs_xc() const { return ((c_functionals_.size() + x_functionals_.size()) > 0); }
    bool needs_vv10() const {return needs_vv10_; };
    bool needs_grac() const {return needs_grac_; };
    bool is_unpolarized() const;
    bool is_meta() const;
    bool is_gga() const;
    bool is_x_lrc() const { return x_omega_ != 0.0; }
    bool is_c_lrc() const { return c_omega_ != 0.0; }
    bool is_x_hybrid() const { return x_alpha_ != 0.0; }
    bool is_c_hybrid() const { return c_alpha_ != 0.0; }
    bool is_c_scs_hybrid() const { return c_os_alpha_ != 0.0 || c_ss_alpha_ != 0.0; }

    // => Utility <= //

    void print(std::string out_fname = "outfile", int print = 1) const;
    void py_print() const { print("outfile", 1); }
    void py_print_detail(int level) const { print("outfile", level); }

};

}

#endif
